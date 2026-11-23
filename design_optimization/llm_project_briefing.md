# LLM Briefing Document: crazieDrop Electropermanent/Flux-Choking Dropper

> **Audience Note (For LLM Consumption):**
> This document is a structured technical briefing designed to give a Large Language Model complete, exact context on the **crazieDrop** hardware, electromechanical baseline, PCB circuit constraints, and optimization problem space without needing to parse raw KiCad files. All values below are verified directly against `crazieDrop.kicad_sch`, `crazieDrop.kicad_pcb`, `drop_test_code/src/main.cpp`, and user specifications.

---

## 1. Mission & System Requirements

- **Platform:** Crazyflie micro-quadrotor (powered from a **1S Crazyflie LiPo / 3.3 V rail**; nominal 3.0–4.2 V under load, ~3.3 V typical under flight load).
- **Payload ("Zippy"):** A **27 g** walking robot ($F_{\text{gravity}} \approx 0.265\text{ N}$, requiring significantly higher passive holding force—typically $\ge 0.6\text{–}1.0\text{ N}$—to survive flight accelerations, prop wash, and off-axis peel/vibration).
- **Total Dropper Mass Budget:** **$< 7.0\text{ g}$** (strict upper bound for the entire dropper mechanism, including PCB/electronics if counted, magnet, steel core, copper coil, and bobbin/housing).
- **Operating Principle (Passive Hold / Active Choke Drop):**
  - **Hold (0 W power):** A permanent magnet sits at the top of a low-carbon steel rod. Flux travels down the steel rod and bridges across into a mating steel target rod on Zippy, holding Zippy passively during flight with zero battery drain.
  - **Drop (Pulsed power):** A copper solenoid wound around the dropper's steel rod is pulsed via a charged capacitor bank so its MMF ($N \cdot I$) opposes the permanent magnet's flux inside the core/air-gap, driving net tip flux $B_{\text{gap}} \to 0$ long enough for gravity to separate Zippy cleanly without re-attracting or reverse-magnetizing the target.

---

## 2. Current Mechanical & Magnetic Baseline

From top (drone side) to bottom (payload side), coaxial cylindrical stack:

1. **PCB Substrate (`crazieDrop.kicad_pcb`):**
   - Dimensions: **18.0 mm × 22.0 mm × 1.6 mm**, 2-layer FR4 (1 oz / 35 µm Cu).
   - Bottom side (`B.Cu` / `B.SilkS`) is completely bare of components and has a **Ø6.0 mm circular alignment silkscreen (`MAG`)** centered at `(149.95, 100.0)`—exact center of the 18 × 22 mm board `(141..159, 89..111)`.
2. **Permanent Magnet (Top of stack, flush against bottom of PCB — Verified McMaster Spec `magnet.png`):**
   - Geometry: **Ø6.35 mm (1/4") × 2.54 mm (1/10") tall** cylinder, Nickel-plated (`±0.004"` tolerances).
   - Grade & Rating: **Grade N52 NdFeB** ($B_r \approx 1.46\text{ T}$, max operating temp **175 °F / 80 °C**), rated **2.8 lbs (12.45 N / 1,270 gf)** direct pull.
   - Estimated mass: **~0.60 g** ($\rho \approx 7.5\text{ g/cm}^3$).
   - **Empirical Bench Observation:** In physical bench testing, the prototype passively held a **60 g rod** (~0.59 N; leftover ~9.5" length of the same 1/4" McMaster magnetic iron rod).
3. **Dropper Core Rod (Stacked directly beneath magnet, coaxial — Verified McMaster Spec `rod.png`):**
   - Material: **Annealed Multipurpose Low-Carbon Magnetic Iron** (consumable-electrode-melted pure magnetic iron, Rockwell B45, **relative magnetic permeability $\mu_r = 4,000$**, $B_{\text{sat}} \approx 2.15\text{ T}$, solid conductive metal $\sigma \approx 1.0 \times 10^7\text{ S/m}$).
   - Geometry: **Ø6.35 mm (1/4") × 9.5 mm long**, careful lathe-cut pole face (~15 µm effective contact roughness air gap).
   - Estimated mass: **~2.37 g** (McMaster 1-ft bar weight = 0.167 lbs = 75.75 g $\implies 0.2485\text{ g/mm}$).
4. **Electromagnet Coil (Wound around the 9.5 mm iron core inside a bobbin/housing):**
   - Current wire gauge: **32 AWG** enameled magnet wire ($d_{\text{bare}} \approx 0.202\text{ mm}$, $d_{\text{ins}} \approx 0.224\text{ mm}$, $\mathcal{R} \approx 0.538\ \Omega/\text{m}$ at 20 °C).
   - Schematic placeholder value: `1 mH` (`L1`), connected to the PCB via 2-pin 2.54 mm through-hole pads.
5. **Target ("Zippy Side" Pickup Slug):**
   - Material: **Same McMaster Annealed Low-Carbon Magnetic Iron ($\mu_r = 4,000$)** (note: an older CAD render showed a second permanent magnet, which is **no longer used** in the baseline—though re-evaluating a second magnet with saturated push-off is an optimization research topic).
   - Geometry: **Ø6.35 mm × 3.0 mm long** (mass ~0.75 g, mounted on Zippy, so this 0.75 g is part of Zippy's 27 g, not the dropper's 7 g budget).

> **Critical Physics Note for Optimization ($F \propto B_{\text{gap}}^2$):**
> Because Zippy's target is **unmagnetized soft low-carbon steel** (not a permanent magnet), magnetic force is **strictly attractive regardless of field polarity**:
> $$F_{\text{mag}} \propto \int_{\text{face}} B_z^2 \, dA$$
> If the coil under-drives ($NI < NI_{\text{cancel}}$), residual PM flux still holds Zippy. If the coil **over-drives** ($NI \gg NI_{\text{cancel}}$), the coil's *reversed* flux magnetizes Zippy's steel slug in the opposite direction and **still attracts it**! Therefore, successful release requires hitting a "null window" around $B_{\text{gap}} \approx 0$ (or sweeping/decaying through zero slowly enough—or holding near zero—for Zippy to fall $x = \frac{1}{2} g t^2$ out of the short-range magnetic capture zone, e.g., $t \approx 10\text{–}20\text{ ms} \implies \Delta x \approx 0.5\text{–}2.0\text{ mm}$).

---

## 3. Complete PCB Circuit Architecture & Analysis

*(Note: `production/bom.csv` reflects an older Revision 1 direct-battery board. The active design in `crazieDrop.kicad_sch` and `crazieDrop.kicad_pcb` is Revision 2, detailed below.)*

### 3.1 Netlist & Component Summary

| Ref | Value / Part | Package | LCSC | Role in Circuit |
|---|---|---|---|---|
| `BT1` | ~3.3 V Crazyflie Rail | `PinHeader_1x02` | — | Input power (`Net-(BT1-+)` and `GND`) |
| `C3` | 22 µF | 0201 (`0603Metric`) | C45783 | Input decoupling cap on `BT1+` (note: 22 µF in 0201 has severe DC-bias/voltage derating; check voltage rating) |
| `U1` | `TLV61046ADB` | SOT-23-6 | — | 28 V max synchronous/diode-integrated step-up DC-DC boost converter (`VI` & `EN` tied to `BT1+`) |
| `L2` | 4.7 µH | Vishay IFSC-1515AH (4×4×1.8 mm) | — | Boost power inductor between `BT1+` and `U1.SW` |
| `R5` | `R_Trim` (Bourns 3314J, 1 MΩ) | 4.5×4.5 mm SMD Pot | C56560 | Top feedback resistor from `V_CAP` (`Net-(D1-K)`) to `U1.FB` |
| `R3` | 30.1 kΩ | 0201 (`0603Metric`) | C2907025 | Bottom feedback resistor from `U1.FB` to `GND` |
| `D3` | `B5819W` (40 V, 1 A Schottky) | SOD-123 | C8598 | Bypass/pre-charge diode from `BT1+` (anode) to `V_CAP` (cathode) |
| `C1` | 1000 µF | Radial THT Ø10 mm, P5.0 mm | — | Bulk energy-storage capacitor on `V_CAP` (`Net-(D1-K)`) to `GND` |
| `L1` | Electromagnet Coil | `PinHeader_1x02` | — | Connected between `V_CAP` (`Net-(D1-K)`) and `Q1` Drain (`Net-(D1-A)`) |
| `D1` | `SS14` (40 V, 1 A Schottky) | SMA (`DO-214AC`) | C2480 | Flyback diode across `L1` (Anode at `Q1.D`, Cathode at `V_CAP`) |
| `Q1` | `AO3400A` (30 V NMOS) | SOT-23 | C20917 | Low-side coil switch (`D` = `L1-2`, `S` = `GND`, `G` = `Net-(D2-K)`) |
| `J1` | `control_gpio` | `PinHeader_1x01` | — | Single-pin logic trigger input from Crazyflie deck GPIO (~3.0–3.3 V logic) |
| `C2` | 2.2 µF | 0603 | C23630 | Series AC-coupling (differentiator) capacitor from `J1` to `R2` |
| `R2` | 100 Ω | 0603 | C22775 | Series gate current-limiting / damping resistor (`C2` to `Q1.G`) |
| `R1` | 100 kΩ | 0603 | C25803 | Gate pull-down / `C2` discharge resistor (`Q1.G` to `GND`) |
| `D2` | `B5819W` (40 V Schottky) | SOD-123 | C8598 | Negative-edge recovery clamp (`Anode` = `GND`, `Cathode` = `Q1.G`) |

---

### 3.2 Subsystem Design Rationale ("Why It Was Built This Way")

#### A. Capacitive Pulsed-Power Architecture (`U1`, `L2`, `R5`, `R3`, `D3`, `C1`)
1. **Decoupling Peak Drop Power from the 3.3 V Crazyflie Battery:**
   - Driving a small, fine-wire electromagnet directly from a 3.3 V micro-drone battery under load causes two problems: (a) 3.3 V limits $V/R$ (and $dI/dt = V/L$) unless wire gauge is very thick, and (b) drawing several amps directly from the Crazyflie 1S LiPo sags the 3.3 V rail and brown-outs the flight controller / motors mid-hover.
   - Instead, `TLV61046ADB` (`U1`) slowly trickles energy from the 3.3 V rail (`I_LIM` ≈ 0.95 A switch current limit) into a large **1000 µF reservoir capacitor (`C1`)** prior to the drop. During the 10–100 ms drop event, practically all coil current is supplied locally by `C1` with minimal transient load on the drone battery.
2. **Tunable Drop Voltage via Trim Pot (`R5` + `R3`) — Crucial for Soft-Steel Nulling:**
   - Because Zippy uses a **passive steel target**, exact flux cancellation requires tuning the coil current $I_{\text{coil}} \approx V_{\text{cap}} / R_{\text{coil}}$ so $N I_{\text{coil}}$ matches the permanent magnet's operating point (neither under-canceling nor over-driving into reverse attraction).
   - `TLV61046A` has $V_{\text{FB}} = 0.795\text{ V}$. With `R3` = 30.1 kΩ and `R5` = 0–1 MΩ trim pot:
     $$V_{\text{out}} = 0.795\text{ V} \times \left(1 + \frac{R_5}{30.1\text{ k}\Omega}\right)$$
     This allows continuous hardware tuning from **$V_{\text{in}}$ (~3.0 V) up to ~27.2 V** (just under `TLV61046A`'s 28 V OVP and `AO3400A`'s 30 V $V_{\text{DS,max}}$) using a screwdriver during drop testing!
3. **Pre-Charge / Bypass Diode (`D3`):**
   - `TLV61046A` includes an internal input-to-output isolation switch, meaning at startup `C1` (1000 µF) is at 0 V. `D3` (`B5819W`) provides a direct path from `BT1+` to `C1` so the cap can pre-charge to ~3.0 V or operate in unboosted direct-battery mode if the boost converter is unpopulated/disabled, and relieves startup stress on `U1`'s internal PMOS isolation FET when charging a massive 1000 µF load. *(Caveat: see Section 3.3 regarding inrush current.)*

#### B. Hardware Fail-Safe One-Shot Gate Driver (`J1`, `C2`, `R2`, `R1`, `D2`, `Q1`)
1. **Why AC-Coupled Gate Drive?**
   - If software crashes, a GPIO defaults HIGH at boot, or a user holds a button too long, continuously connecting a low-resistance coil across `C1` / `U1` would burn up the 32 AWG coil, melt the 3D-printed bobbin, and continuously overload the boost converter / battery.
   - Series capacitor `C2` (2.2 µF) + pull-down `R1` (100 kΩ) form a **hardware high-pass differentiator (monostable pulse limiter)**:
     - On a rising edge at `J1` ($0 \to 3.3\text{ V}$), `C2` instantaneously steps `Q1`'s gate to ~3.3 V (well above `AO3400A`'s $V_{\text{GS(th)}} \approx 1.05\text{ V}$, giving $R_{\text{DS(on)}} \approx 28\text{–}35\text{ m}\Omega$).
     - Even if `J1` stays stuck HIGH forever, `C2` charges through `R1` with time constant:
       $$\tau_{\text{gate}} = R_1 C_2 = 100\text{ k}\Omega \times 2.2\ \mu\text{F} = 220\text{ ms}$$
       Gate voltage decays as $V_G(t) = 3.3\,e^{-t / 0.22\text{ s}}$, automatically shutting off `Q1` after ~200–250 ms ($V_G \to V_{\text{th}}$).
     - Meanwhile, firmware can still command *shorter* pulses (e.g., 10–100 ms) simply by pulling `J1` LOW early.
2. **Why `D2` (`B5819W`) Is on the Gate:**
   - When `J1` goes LOW ($3.3\text{ V} \to 0\text{ V}$) after `C2` has partially charged, `Q1`'s gate would otherwise swing negative to $-3.3\text{ V}$ and take ~500 ms to recover through `R1`, preventing rapid re-triggering. Schottky `D2` clamps the negative swing at $-0.3\text{ V}$ and rapidly resets `C2` via `R2` ($\tau_{\text{reset}} = R_2 C_2 = 100\ \Omega \times 2.2\ \mu\text{F} = 0.22\text{ ms}$) so the circuit is immediately ready for another pulse.

#### C. Physical Layout & Mechanical Integration (`crazieDrop.kicad_pcb`)
- **All components are single-sided on `F.Cu` (Top Layer):** Keeps `B.Cu` completely flat so the Ø6.35 mm permanent magnet (`MAG` silkscreen at board center) can glue/sit flush against the underside of the PCB, making the PCB itself the structural top cap of the dropper assembly (as shown in the CAD render).
- **Trace Sizing & Grounding:** Power paths (`BT1+`, `V_CAP`, `Q1.D`, `GND`) use 0.8–1.0 mm wide traces on `F.Cu` backed by a full `B.Cu` ground pour.

---

### 3.3 Electrical Quirks & Constraints Relevant to Optimization

1. **Capacitor `C1` Voltage Rating vs. Mass/Size Tradeoff:**
   - A 1000 µF radial electrolytic in Ø10 × 16 mm (`CP_Radial_D10.0mm_P5.00mm`) is typically rated for **6.3 V, 10 V, or 16 V** (25 V / 35 V 1000 µF caps are usually Ø10×20 mm or Ø12.5×20 mm and weigh **2.0–3.5 g alone!**).
   - Therefore, the actual usable $V_{\text{cap}}$ is bounded by `C1`'s voltage rating ($V_{\text{rated}} \le 28\text{ V}$) and `C1`'s physical mass within the 7 g budget. If optimization reduces required pulse energy $E_{\text{pulse}}$, shrinking `C1` saves substantial mass.
2. **Capacitor Discharge Dynamics (RLC / RC Decay, Not Constant Voltage!):**
   - Because `U1` (`TLV61046A`) can only supply ~50–150 mA at boosted output voltages while the coil pulse draws ~0.5–5 A, **`C1` acts as a finite energy tank** during the drop pulse:
     $$E_{\text{stored}} = \frac{1}{2} C_1 V_0^2$$
   - The coil + capacitor + MOSFET form an **overdamped/underdamped RLC discharge** (clamped from oscillating negative by `D1` & `D3`):
     $$L_{\text{coil}} \frac{dI}{dt} + (R_{\text{coil}} + R_{\text{ESR}} + R_{\text{DS(on)}}) I = V_C(t), \qquad C_1 \frac{dV_C}{dt} = -I(t) + I_{\text{boost}}$$
   - This is actually **advantageous** if tuned properly: $I(t)$ rises to peak near $I_{\text{cancel}}$ and naturally droops as `C1` discharges, or can be sized so `C1` holds $V_C$ nearly flat over a 15–30 ms release window.
3. **Slow Flyback Decay (`D1` = `SS14`):**
   - Standard freewheeling diode `D1` clamps turn-off voltage to ~0.4 V, meaning coil current decays with $\tau = L / R_{\text{coil}}$ after `Q1` opens. That keeps opposing flux active slightly longer after GPIO turn-off without high-voltage inductive kickback on the 30 V `AO3400A`.
4. **Solid Steel Core Eddy Currents:**
   - Because the Ø6.35 mm low-carbon steel rod is solid metal (not laminated or ferrite), rapid $dB/dt$ from a fast coil pulse induces circumferential eddy currents in the rod that **shield the interior of the core** and delay flux cancellation at the pole face!
   - Magnetic diffusion time constant for a cylinder of radius $a = 3.175\text{ mm}$:
     $$\tau_{\text{eddy}} \approx \frac{\mu_0 \mu_r \sigma a^2}{(2.405)^2}$$
     For low-carbon steel ($\sigma \approx 6\times 10^6\text{ S/m}$, $\mu_r \approx 100\text{–}300$), $\tau_{\text{eddy}} \approx 1.3\text{–}4.0\text{ ms}$ (and full penetration takes ~3$\tau \approx 4\text{–}12\text{ ms}$). Very short sub-millisecond high-voltage pulses will waste energy in eddy heating without choking the central core flux; pulse widths on the order of **10–40 ms** (and/or thinner/hollow/slotted core geometries) penetrate much more effectively.
5. **Minor PCB Notes (Do Not Touch KiCad Now, Just Context):**
   - `TLV61046A` datasheet specifies $L = 2.2\text{–}10\ \mu\text{H}$ and output capacitance typically $1\text{–}22\ \mu\text{F}$ ceramic; driving a 1000 µF electrolytic directly on `VOUT` can cause `U1` to hit cycle-by-cycle current limit / startup fault timeout if `U1` tries to boost from 3.3 V to high voltage all at once, though `D3` pre-charges `C1` to ~3.0 V at plug-in (watch battery connector hot-plug spark/inrush through `D3`).
   - If `J1` is left held HIGH and `Q1` turns off via `R1*C2` decay (~220 ms), `Q1`'s gate passes slowly through its linear region ($V_{\text{GS}} \approx 1.0\text{–}1.8\text{ V}$); driving `J1` with a clean digital pulse (e.g., 15–50 ms in firmware) avoids linear-region FET heating while keeping `R1*C2` purely as a fault safety net.

---

## 4. Optimization Problem Formulation Summary

When optimizing the **electromagnet dropper** for the existing PCB:

### Fixed / Hard Constraints
- **Battery input:** $V_{\text{bat}} \approx 3.3\text{ V}$ (Crazyflie 1S LiPo under hover load).
- **Drive voltage range:** $V_0 \in [3.0\text{ V},\ 27\text{ V}]$ (set via `R5` trim pot; bounded by `C1` voltage rating and `U1` 28 V / `Q1` 30 V limits).
- **Bulk capacitance:** $C_1 = 1000\ \mu\text{F}$ nominal on current PCB (or smaller THT/SMD equivalent fitting 5 mm THT pads if mass-optimized).
- **Mass limit:** $m_{\text{total}} < 7.0\text{ g}$ (track magnet + core + copper + bobbin/housing + PCB/cap breakdown explicitly!).
- **Holding requirement:** Reliably hold **27 g Zippy** ($F_{\text{hold}} \ge SF \times 0.265\text{ N}$, where safety factor $SF \approx 2.5\text{–}4\times \implies F_{\text{hold}} \approx 0.7\text{–}1.1\text{ N}$).
- **Release requirement:** During coil energization, reduce attractive force at the dropper-to-Zippy interface below $m_{\text{Zippy}} g = 0.265\text{ N}$ (ideally $< 0.05\text{ N}$) for long enough ($\Delta t_{\text{null}} \ge 15\text{–}25\text{ ms}$) that Zippy falls clear of the magnetic capture range **without over-driving $B_{\text{gap}}$ negative and re-attracting the soft-steel slug**.

### Primary Design Variables to Optimize
1. **Permanent Magnet Geometry & Grade:** Diameter $D_m$, thickness $h_m$, grade (N42/N52), and axial position/flux-shunt geometry.
2. **Dropper Core Geometry:** Diameter $D_c$, length $L_c$, taper/pole-tip geometry, or hollow/stepped cross-section (trading off passive flux guiding vs. mass $\propto D_c^2 L_c$ vs. eddy-current diffusion time $\propto D_c^2$).
3. **Target ("Zippy") Slug Geometry:** Diameter $D_t$ and thickness $h_t$ (currently Ø6.35 × 3.0 mm low-carbon steel; thinner slugs saturate earlier, which caps holding force).
4. **Coil Winding & Wire Gauge:** AWG gauge (e.g., 28–38 AWG), number of turns $N$, axial winding length & radial build, coil resistance $R_{\text{coil}}$, inductance $L_{\text{coil}}$, and axial placement along the rod (near tip vs. full length—winding closer to the pole tip chokes gap flux with less leakage and less demagnetizing H-field on the permanent magnet!).
5. **Electrical Operating Point:** Boost target voltage $V_0$ (via `R5` trim pot) and GPIO pulse width $t_{\text{pulse}}$ so $I(t)$ hits the null window cleanly within `C1`'s stored energy.
