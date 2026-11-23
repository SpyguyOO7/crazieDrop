from dataclasses import dataclass
import numpy as np
from shapely.geometry import Point, Polygon
from netgen.occ import OCCGeometry, gp_Pnt, Segment, Wire, Face, Glue
import ngsolve as ngs

from materials import Material, WireSpec, AIR, CU_DENSITY_KG_M3


# --- 2D Axisymmetric Polygon Part ---
# Coordinates are strictly (r_mm, z_mm) in millimeters with r_mm >= 0.

@dataclass(frozen=True)
class AxiPart:
    name: str
    pts_rz_mm: tuple[tuple[float, float], ...]
    mat: Material
    maxh_mm: float = 0.4
    group: str = "dropper"  # "dropper", "target", or "static"
    wire: WireSpec | None = None
    turns: float | None = None
    coil_sign: float = 1.0  # +1.0 = +phi current (+z field), -1.0 = -phi current (-z opposing field)

    def __post_init__(self):
        if len(self.pts_rz_mm) < 3:
            raise ValueError(f"Part '{self.name}' needs at least 3 (r, z) points.")
        if self.group not in ("dropper", "target", "static"):
            raise ValueError(f"Part '{self.name}' group must be 'dropper', 'target', or 'static'.")
        if self.coil_sign not in (-1.0, 1.0):
            raise ValueError(f"Part '{self.name}' coil_sign must be +1.0 or -1.0.")

        pts = np.asarray(self.pts_rz_mm, dtype=float)
        if pts.ndim != 2 or pts.shape[1] != 2:
            raise ValueError(f"Part '{self.name}' points must have shape (N, 2).")
        if np.any(pts[:, 0] < -1e-9):
            raise ValueError(f"Part '{self.name}' has negative radius r < 0: min r = {pts[:, 0].min()}.")

        poly = Polygon(self.pts_rz_mm)
        if not poly.is_valid or poly.area <= 1e-9:
            raise ValueError(f"Part '{self.name}' forms an invalid or zero-area polygon.")

        if self.mat.kind == "coil":
            if self.wire is None and self.turns is None:
                raise ValueError(f"Coil part '{self.name}' requires either 'wire' or 'turns'.")

    @property
    def poly_mm(self) -> Polygon:
        return Polygon(self.pts_rz_mm)

    @property
    def area_mm2(self) -> float:
        return float(self.poly_mm.area)

    @property
    def centroid_r_mm(self) -> float:
        return float(self.poly_mm.centroid.x)

    @property
    def volume_mm3(self) -> float:
        # Exact volume of revolution by Pappus's centroid theorem: V = 2 * pi * r_centroid * Area
        return 2.0 * np.pi * self.centroid_r_mm * self.area_mm2

    @property
    def n_turns(self) -> float:
        if self.mat.kind != "coil":
            return 0.0
        if self.turns is not None:
            return float(self.turns)
        assert self.wire is not None
        # Packing definition: N * area_bare_cu = fill_factor * window_area
        return (self.wire.fill_factor * self.area_mm2) / self.wire.area_cu_mm2

    @property
    def coil_res_ohm(self) -> float:
        if self.mat.kind != "coil" or self.wire is None:
            return 0.0
        mean_turn_len_m = 2.0 * np.pi * (self.centroid_r_mm * 1e-3)
        total_len_m = self.n_turns * mean_turn_len_m
        return total_len_m * self.wire.ohm_per_m

    @property
    def mass_g(self) -> float:
        vol_m3 = self.volume_mm3 * 1e-9
        if self.mat.kind == "coil":
            if self.wire is not None:
                mean_turn_len_m = 2.0 * np.pi * (self.centroid_r_mm * 1e-3)
                cu_vol_m3 = self.n_turns * mean_turn_len_m * (self.wire.area_cu_mm2 * 1e-6)
                return cu_vol_m3 * CU_DENSITY_KG_M3 * 1e3
            return 0.0
        return vol_m3 * self.mat.density_kg_m3 * 1e3

    def with_z_shift(self, dz_mm: float) -> "AxiPart":
        shifted_pts = tuple((r, z + dz_mm) for r, z in self.pts_rz_mm)
        return AxiPart(
            name=self.name,
            pts_rz_mm=shifted_pts,
            mat=self.mat,
            maxh_mm=self.maxh_mm,
            group=self.group,
            wire=self.wire,
            turns=self.turns,
            coil_sign=self.coil_sign,
        )


def rect_rz(
    name: str,
    r_in_mm: float,
    r_out_mm: float,
    z_bot_mm: float,
    z_top_mm: float,
    mat: Material,
    maxh_mm: float = 0.4,
    group: str = "dropper",
    wire: WireSpec | None = None,
    turns: float | None = None,
    coil_sign: float = 1.0,
) -> AxiPart:
    if r_out_mm <= r_in_mm:
        raise ValueError(f"Part '{name}': r_out_mm ({r_out_mm}) must exceed r_in_mm ({r_in_mm}).")
    if z_top_mm <= z_bot_mm:
        raise ValueError(f"Part '{name}': z_top_mm ({z_top_mm}) must exceed z_bot_mm ({z_bot_mm}).")

    pts = (
        (r_in_mm, z_bot_mm),
        (r_out_mm, z_bot_mm),
        (r_out_mm, z_top_mm),
        (r_in_mm, z_top_mm),
    )
    return AxiPart(
        name=name,
        pts_rz_mm=pts,
        mat=mat,
        maxh_mm=maxh_mm,
        group=group,
        wire=wire,
        turns=turns,
        coil_sign=coil_sign,
    )


COIL_MAT = Material(name="copper_coil", kind="coil")


# --- Assembly Validation & OCC Mesh Generation ---

@dataclass
class AxiModel:
    parts: list[AxiPart]
    mesh: ngs.Mesh
    air_r_max_mm: float
    air_z_min_mm: float
    air_z_max_mm: float
    mst_box_mm: tuple[float, float, float, float] | None  # (0.0, r_box, z_bot, z_top)
    mst_shell_mm: tuple[float, float, float, float, float, float] | None = None  # (r_in, r_out, z_t_min, z_bot, z_t_max, z_top)

    def mass_breakdown_g(self, group: str = "dropper") -> dict[str, float]:
        out = {p.name: p.mass_g for p in self.parts if p.group == group}
        out["total"] = sum(out.values())
        return out


def validate_parts(parts: list[AxiPart]) -> None:
    names = [p.name for p in parts]
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate part names detected: {names}")

    # Pairwise overlap check (shared edges/vertices are allowed; overlapping interiors are forbidden)
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            inter_area = parts[i].poly_mm.intersection(parts[j].poly_mm).area
            if inter_area > 1e-7:
                raise ValueError(
                    f"Geometry overlap error: '{parts[i].name}' and '{parts[j].name}' "
                    f"intersect by {inter_area:.5f} mm^2."
                )


def _poly_to_occ_face_m(pts_rz_mm: tuple[tuple[float, float], ...], name: str, maxh_mm: float) -> Face:
    # Ensure counter-clockwise vertex ordering and convert mm -> meters at mesh construction time
    poly = Polygon(pts_rz_mm)
    coords = list(poly.exterior.coords)[:-1]
    if not poly.exterior.is_ccw:
        coords = list(reversed(coords))

    pnts = [gp_Pnt(max(0.0, r) * 1e-3, z * 1e-3, 0.0) for r, z in coords]
    segs = [Segment(pnts[k], pnts[(k + 1) % len(pnts)]) for k in range(len(pnts))]
    face = Face(Wire(segs))
    face.faces.name = name
    face.faces.maxh = maxh_mm * 1e-3

    # Label symmetry axis edges (r == 0) vs interior edges
    for edge in face.edges:
        c_r = edge.center.x
        if abs(c_r) < 1e-9:
            edge.name = "axis"
        else:
            edge.name = f"if_{name}"
    return face


def build_axi_model(
    parts: list[AxiPart],
    air_pad_r_mm: float = 25.0,
    air_pad_z_mm: float = 25.0,
    air_maxh_mm: float = 2.5,
    gap_refine_mm: float = 0.12,
) -> AxiModel:
    validate_parts(parts)

    all_r = [r for p in parts for (r, _) in p.pts_rz_mm]
    all_z = [z for p in parts for (_, z) in p.pts_rz_mm]
    r_max_mm = max(all_r) + air_pad_r_mm
    z_min_mm = min(all_z) - air_pad_z_mm
    z_max_mm = max(all_z) + air_pad_z_mm

    # Optional dedicated Maxwell Stress Tensor air box enclosing all "target" parts
    target_parts = [p for p in parts if p.group == "target"]
    mst_box_mm = None
    mst_shell_mm = None

    part_faces = [_poly_to_occ_face_m(p.pts_rz_mm, p.name, p.maxh_mm) for p in parts]

    # Air domain (in meters)
    air_pts_mm = ((0.0, z_min_mm), (r_max_mm, z_min_mm), (r_max_mm, z_max_mm), (0.0, z_max_mm))
    air_face = _poly_to_occ_face_m(air_pts_mm, "air_outer", air_maxh_mm)
    for edge in air_face.edges:
        if abs(edge.center.x) < 1e-9:
            edge.name = "axis"
        else:
            edge.name = "outer"

    # If target parts exist, build a clean mid-gap air buffer box around target for accurate MST integration
    extra_air_faces = []
    gap_mm_val = None
    gap_z_mid_mm = None
    if target_parts:
        t_r_max = max(r for p in target_parts for (r, _) in p.pts_rz_mm)
        t_z_min = min(z for p in target_parts for (_, z) in p.pts_rz_mm)
        t_z_max = max(z for p in target_parts for (_, z) in p.pts_rz_mm)

        non_target = [p for p in parts if p.group != "target"]
        if non_target:
            nt_z_min = min(z for p in non_target for (_, z) in p.pts_rz_mm if z >= t_z_max - 1e-6) if any(
                z >= t_z_max - 1e-6 for p in non_target for (_, z) in p.pts_rz_mm
            ) else t_z_max + 2.0
            gap_mm = nt_z_min - t_z_max
            if gap_mm >= 0.004:
                gap_mm_val = gap_mm
                z_top_box = t_z_max + 0.80 * gap_mm
                gap_z_mid_mm = t_z_max + 0.50 * gap_mm
                r_box = t_r_max + max(1.5, 0.5 * gap_mm)
                z_bot_box = t_z_min - max(1.5, 0.5 * gap_mm)
                mst_box_mm = (0.0, r_box, z_bot_box, z_top_box)
                mst_shell_mm = (t_r_max, r_box, t_z_min, z_bot_box, t_z_max, z_top_box)

                box_maxh = min(0.35, max(gap_refine_mm, 0.4 * gap_mm))
                box_pts = ((0.0, z_bot_box), (r_box, z_bot_box), (r_box, z_top_box), (0.0, z_top_box))
                mst_face = _poly_to_occ_face_m(box_pts, "air_mst", box_maxh)
                for pf in part_faces:
                    mst_face = mst_face - pf
                mst_face.faces.name = "air_mst"
                mst_face.faces.maxh = box_maxh * 1e-3
                extra_air_faces.append(mst_face)

    # Subtract all parts and inner air boxes from outer air
    air_rem = air_face
    for pf in part_faces:
        air_rem = air_rem - pf
    for ef in extra_air_faces:
        air_rem = air_rem - ef
    air_rem.faces.name = "air_outer"
    air_rem.faces.maxh = air_maxh_mm * 1e-3

    glued = Glue([air_rem] + extra_air_faces + part_faces)

    # Re-assign face names AFTER Boolean cuts + Glue by matching both centroid AND area
    # (prevents non-convex air faces with holes whose centroid falls inside a hole from matching the hole part!)
    for face in glued.faces:
        fc_r_mm = face.center.x * 1e3
        fc_z_mm = face.center.y * 1e3
        f_area_mm2 = face.mass * 1e6
        p_test = Point(fc_r_mm, fc_z_mm)

        matched_part = None
        for p in parts:
            if p.poly_mm.contains(p_test) and abs(f_area_mm2 - p.area_mm2) < 0.02 * p.area_mm2 + 1e-4:
                matched_part = p
                break

        if matched_part is not None:
            face.name = matched_part.name
            face.maxh = matched_part.maxh_mm * 1e-3
        elif mst_box_mm is not None and f_area_mm2 < (mst_box_mm[1] * (mst_box_mm[3] - mst_box_mm[2]) + 1e-3) and (
            mst_box_mm[0] <= fc_r_mm <= mst_box_mm[1]
            and mst_box_mm[2] <= fc_z_mm <= mst_box_mm[3]
        ):
            face.name = "air_mst"
        else:
            face.name = "air_outer"

    # Label boundaries AFTER Glue so split/shared edges never lose or overwrite names
    mst_r = mst_box_mm[1] if mst_box_mm is not None else -1.0
    mst_z_bot = mst_box_mm[2] if mst_box_mm is not None else 0.0
    mst_z_top = mst_box_mm[3] if mst_box_mm is not None else 0.0

    for edge in glued.edges:
        c_r = edge.center.x * 1e3
        c_z = edge.center.y * 1e3
        if abs(c_r) < 1e-5:
            edge.name = "axis"
            if z_min_mm + 2.0 < c_z < z_max_mm - 2.0:
                edge.maxh = min(0.35, air_maxh_mm) * 1e-3
        elif (
            abs(c_r - r_max_mm) < 1e-4
            or abs(c_z - z_min_mm) < 1e-4
            or abs(c_z - z_max_mm) < 1e-4
        ):
            edge.name = "outer"
        elif mst_box_mm is not None:
            if abs(c_z - mst_z_top) < 1e-4 and c_r <= mst_r + 1e-4:
                edge.name = "mst_top"
                if gap_mm_val is not None:
                    edge.maxh = min(0.12, max(0.025, 1.5 * gap_mm_val)) * 1e-3
            elif abs(c_z - mst_z_bot) < 1e-4 and c_r <= mst_r + 1e-4:
                edge.name = "mst_bot"
            elif abs(c_r - mst_r) < 1e-4 and mst_z_bot - 1e-4 <= c_z <= mst_z_top + 1e-4:
                edge.name = "mst_side"

        # Also refine the pole-face and target-face edges directly adjacent to tiny gaps
        if gap_mm_val is not None and gap_z_mid_mm is not None:
            if abs(c_z - gap_z_mid_mm) < (0.6 * gap_mm_val + 1e-3) and c_r <= mst_r + 1e-4:
                edge.maxh = min(0.12, max(0.025, 1.5 * gap_mm_val)) * 1e-3

    geo = OCCGeometry(glued, dim=2)
    ngmesh = geo.GenerateMesh()
    mesh = ngs.Mesh(ngmesh)

    return AxiModel(
        parts=parts,
        mesh=mesh,
        air_r_max_mm=r_max_mm,
        air_z_min_mm=z_min_mm,
        mst_shell_mm=mst_shell_mm,
        air_z_max_mm=z_max_mm,
        mst_box_mm=mst_box_mm,
    )
