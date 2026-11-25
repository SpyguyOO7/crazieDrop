# Schematic Review & JLCPCB / LCSC Component Findings

> Verified against live JLCPCB / LCSC API (`scripts/query_jlc_parts.py`) on 2026-11-25, aligned to the current `crazieDrop.kicad_sch` reference designators (`R6` = reservoir charge resistor, `R4` = series feedback stop, `R3` = bottom feedback resistor).

---

## 1. Key Schematic & BOM Findings

1. **Boost Output Decoupling & Reservoir Isolation (`U1` = `TLV61046ADBVR`):**
   - `TLV61046A` has a **95 ms soft-start / fault latch-off timer** and switches at ~1.05 MHz. Driving `C1` (`1000–1500 µF`) directly on `U1.VOUT` lacks high-frequency MLCC bypassing and risks tripping the 95 ms fault timeout at startup or after a drop pulse.
   - **Fix:** Place local ceramic bypass **`C4` (`10 µF 50 V X5R`, `0805`, `C440198`)** directly from `U1.VOUT` to `GND`, keep the `R5 + R4 + R3` feedback divider on `U1.VOUT`, and isolate the `VCAP` (`C1`) reservoir rail through series charge resistor **`R6` (`330 Ω`, `0805`, `C17630`)**.
   - *(Note on `R6` footprint in schematic:* Make sure `R6` uses `Resistor_SMD:R_0805_2012Metric...` rather than `Capacitor_SMD:C_0805_2012Metric...`.)*
   - At $V_{\text{OUT}} = 20\text{ V}$ and $V_{\text{CAP}} = 3\text{ V}$, `330 Ω` limits peak load on `U1.VOUT` to ~51 mA (safely below `TLV61046A`'s ~80 mA output capability at 20 V from 3.3 V), preventing fault shutdown while recharging `1000 µF` in ~1.0 s (`470 µF` in ~0.47 s, `1500 µF` in ~1.5 s).

2. **3.3 V Battery Tail Handoff (`D3`):**
   - Keep Schottky diode **`D3` (`B5819W`, `SOD-123`, `C8598`)** from `BT1+` (Anode) to `VCAP` (Cathode). Once `C1` discharges down to $V_{\text{BAT}} - V_f \approx 3.0\text{ V}$ (~11 ms into the pulse), `D3` forward-biases and hands off coil tail current (~0.48 A) directly to the 1S LiPo.
   - **BOM Fix:** Change `D3`'s `LCSC` field from `C2480` (`SS14` in `SMA`) to **`C8598`** (`B5819W` in `SOD-123`).

3. **Trim Pot `R5` (`C56560`) Value & Small Feedback Divider Resistors (`R4` + `R3` — `0603` or `0402`):**
   - Live LCSC lookup showed **`C56560` is `3314J-1-504E` (`500 kΩ`)**, **not** `1 MΩ` (`1 MΩ` is `C55491` = `3314J-1-105E`).
   - **Why `R4` & `R3` can be tiny (`0603` or `0402`, NOT `0805`):** The feedback divider only carries $I_{\text{FB}} = 0.795\text{ V} / 24\text{ k}\Omega = \mathbf{33\ \mu\text{A}}$!
     - Worst-case power in **`R4` (`220 kΩ`)** (when `R5 = 0 Ω` at 25 V) is $(24.2\text{ V})^2 / 220\text{ k}\Omega = \mathbf{2.66\text{ mW}}$, and normal steady-state power is **$0.24\text{ mW}$** (with `< 19 V` across `R4`, well below `0402`'s `50 V / 62.5 mW` and `0603`'s `75 V / 100 mW` ratings).
     - Power in **`R3` (`24 kΩ` / `22 kΩ`)** is fixed at $(0.795\text{ V})^2 / 22\text{ k}\Omega = \mathbf{0.029\text{ mW}}$ (`0.795 V` across it).
     - So **`0603`** (`R_0603_1608Metric`, easiest hand-soldering, matches `R1`/`R2`) or **`0402`** (`R_0402_1005Metric`, ultra-compact PCB footprint) both have >20× safety margin—only **`R6` (`330 Ω` charge resistor)** needs `0805`!
   - **Recommended divider with `C56560` (`500 kΩ` pot):**
     - **`R4 = 220 kΩ`** (`0603` `C22961` or `0402` `C25767`) in series between `R5` and `FB`, and **`R3 = 24 kΩ`** (`0603` `C23352` or `0402` `C25769`) from `FB` to `GND` $\implies$ sweeps **`8.08 V – 24.64 V`**.
     - *(Alternative `R3 = 22 kΩ`: `0603` `C31850` or `0402` `C25768` $\implies$ sweeps **`8.75 V – 26.81 V`**.)*
     - *(Alternative if using `1 MΩ` pot `C55491`: `R4 = 330 kΩ` `C23137` + `R3 = 43 kΩ` `C23172` $\implies$ sweeps **`6.89 V – 25.39 V`**.)*
   - Tie unused potentiometer Pad 3 to Wiper Pad 2 on `Bourns 3314J` so wiper bounce under vibration cannot open-circuit `FB`.

4. **KiCad Imperial vs. Metric Footprint Reference:**
   - `0201_0603Metric` = **Imperial `0201`** (`0.6 × 0.3 mm`)
   - `0402_1005Metric` = **Imperial `0402`** (`1.0 × 0.5 mm`)
   - `0603_1608Metric` = **Imperial `0603`** (`1.6 × 0.8 mm`)
   - `0805_2012Metric` = **Imperial `0805`** (`2.0 × 1.2 mm`)

5. **Minor Schematic Hygiene:**
   - Use polarized symbol `Device:C_Polarized` for `C1` (`+` = Pad 1 on `VCAP`, `-` = Pad 2 on `GND`) so polarity is explicit.
   - Add a `PWR_FLAG` on `U1` Pin 1 (`SW`) net (or set Pin 1 passive) to clear KiCad ERC `[power_pin_not_driven]`.

---

## 2. Recommended Complete BOM (Aligned to Current Schematic Designators)

| Ref | Role | Schematic Value | KiCad Footprint (Imperial) | LCSC # | MPN | Brand | JLC Tier | Stock | Price (1+) | Verified Specs |
|---|---|---|---|---|---|---|---|---:|---:|---|
| **`U1`** | Boost IC | `TLV61046ADBVR` | `Package_TO_SOT_SMD:SOT-23-6` | **`C464815`** | `TLV61046ADBVR` | TI | Extended | 12,997 | $0.7102 | 28V Boost, Load Disconnect, SOT-23-6 |
| **`L2`** | Boost Inductor | `4.7uH` | `Inductor_SMD:L_Vishay_IFSC-1515AH_4x4x1.8mm` | **`C354574`** | `CKCS4018-4.7uH/M` | CENKER | Extended | 7,076 | $0.0416 | 4.7µH, 1.7A Isat / 1.2A Irms, 90mΩ, 4×4×1.8mm |
| **`C3`** *(0805)* | Input Bypass | `22uF 25V` | `Capacitor_SMD:C_0805_2012Metric` (**`0805`**) | **`C45783`** | `CL21A226MAQNNNE` | Samsung | **Basic** | 4,738,662 | $0.1084 | 22µF 25V X5R ±20% `0805` *(currently in sch)* |
| **`C3`** *(opt 0603)* | Input Bypass | `10uF 25V` | `Capacitor_SMD:C_0603_1608Metric` (**`0603`**) | **`C96446`** | `CL10A106MA8NRNC` | Samsung | **Basic** | 5,342,707 | $0.0331 | 10µF 25V X5R ±20% `0603` |
| **`C4`** | Local `VOUT` Cap | `10uF 50V` | `Capacitor_SMD:C_0805_2012Metric` (**`0805`**) | **`C440198`** | `GRM21BR61H106KE43L` | Murata | **Basic** | 922,412 | $0.1134 | 10µF **50V** X5R ±10% `0805` |
| **`R6`** | `VCAP` Charge Res | `330R` | `Resistor_SMD:R_0805_2012Metric` (**`0805`**) | **`C17630`** | `0805W8F3300T5E` | UNI-ROYAL | **Basic** | 1,151,059 | $0.0044 | 330Ω 1% 125mW 150V `0805` *(~1.0s charge on 1000µF)* |
| **`R5`** | Trim Pot | `500k` | `Potentiometer_Bourns_3314J_Vertical` | **`C56560`** | `3314J-1-504E` | Bourns | Extended | 3,630 | $0.5924 | **500kΩ** ±20% 250mW 4.5×4.5mm SMD Trim Pot |
| **`R4`** *(0603)* | FB Series Stop | `220k` | `Resistor_SMD:R_0603_1608Metric` (**`0603`**) | **`C22961`** | `0603WAF2203T5E` | UNI-ROYAL | **Basic** | 870,991 | $0.0015 | 220kΩ 1% 100mW 75V `0603` |
| **`R4`** *(opt 0402)* | FB Series Stop | `220k` | `Resistor_SMD:R_0402_1005Metric` (**`0402`**) | **`C25767`** | `0402WGF2203TCE` | UNI-ROYAL | Preferred | 338,183 | $0.0013 | 220kΩ 1% 62.5mW 50V `0402` |
| **`R3`** *(0603, 8.1–24.6V)* | Bottom FB | `24k` | `Resistor_SMD:R_0603_1608Metric` (**`0603`**) | **`C23352`** | `0603WAF2402T5E` | UNI-ROYAL | **Basic** | 632,669 | $0.0019 | 24kΩ 1% 100mW 75V `0603` |
| **`R3`** *(0402, 8.1–24.6V)* | Bottom FB | `24k` | `Resistor_SMD:R_0402_1005Metric` (**`0402`**) | **`C25769`** | `0402WGF2402TCE` | UNI-ROYAL | Preferred | 1,489,995 | $0.0012 | 24kΩ 1% 62.5mW 50V `0402` |
| **`R3`** *(0603, 8.8–26.8V)* | Bottom FB | `22k` | `Resistor_SMD:R_0603_1608Metric` (**`0603`**) | **`C31850`** | `0603WAF2202T5E` | UNI-ROYAL | **Basic** | 2,097,947 | $0.0042 | 22kΩ 1% 100mW 75V `0603` |
| **`R3`** *(0402, 8.8–26.8V)* | Bottom FB | `22k` | `Resistor_SMD:R_0402_1005Metric` (**`0402`**) | **`C25768`** | `0402WGF2202TCE` | UNI-ROYAL | **Basic** | 1,404,966 | $0.0010 | 22kΩ 1% 62.5mW 50V `0402` |
| **`D1`** | Coil Flyback | `SS14` | `Diode_SMD:D_SMA_Handsoldering` | **`C2480`** | `SS14` | MDD | **Basic** | 1,153,154 | $0.0192 | 40V 1A (25A surge) Schottky `SMA` |
| **`D2, D3`** | Gate Clamp / 3.3V Tail | `B5819W` | `Diode_SMD:D_SOD-123` | **`C8598`** | `B5819W SL` | CJ | **Basic** | 951,715 | $0.0248 | 40V 1A (9A surge) Schottky `SOD-123` |
| **`Q1`** | Low-Side NMOS | `AO3400A` | `Package_TO_SOT_SMD:SOT-23_Handsoldering` | **`C20917`** | `AO3400A` | AOS | **Basic** | 1,462,022 | $0.0786 | 30V 5.7A NMOS, $R_{\text{DS(on)}} \le 48\text{ m}\Omega @ 2.5\text{V}$, `SOT-23` |
| **`C2`** | Gate AC-Couple | `2.2uF 16V` | `Capacitor_SMD:C_0603_1608Metric` (**`0603`**) | **`C23630`** | `CL10A225KO8NNNC` | Samsung | **Basic** | 2,504,575 | $0.0097 | 2.2µF 16V X5R ±10% `0603` |
| **`R1`** | Gate Pull-Down | `100k` | `Resistor_SMD:R_0603_1608Metric` (**`0603`**) | **`C25803`** | `0603WAF1003T5E` | UNI-ROYAL | **Basic** | 12,997,360 | $0.0026 | 100kΩ 1% 100mW `0603` |
| **`R2`** | Gate Series | `100R` | `Resistor_SMD:R_0603_1608Metric` (**`0603`**) | **`C22775`** | `0603WAF1000T5E` | UNI-ROYAL | **Basic** | 12,826,755 | $0.0017 | 100Ω 1% 100mW `0603` |

### Alternative `1 MΩ` Potentiometer Divider Set (If Using `1 MΩ` Pot Instead of `C56560`)

| Ref | Schematic Value | KiCad Footprint (Imperial) | LCSC # | MPN | Brand | JLC Tier | Stock | Price (1+) | Notes |
|---|---|---|---|---|---|---|---:|---:|---|
| **`R5`** | `1M` | `Potentiometer_Bourns_3314J_Vertical` | **`C55491`** | `3314J-1-105E` | Bourns | Extended | 841 | $0.6399 | 1MΩ ±20% 250mW 4.5×4.5mm |
| **`R4`** | `330k` | `Resistor_SMD:R_0603_1608Metric` (**`0603`**) | **`C23137`** | `0603WAF3303T5E` | UNI-ROYAL | **Basic** | 742,459 | $0.0017 | Series min-stop |
| **`R3`** | `43k` | `Resistor_SMD:R_0603_1608Metric` (**`0603`**) | **`C23172`** | `0603WAF4302T5E` | UNI-ROYAL | Preferred | 213,410 | $0.0017 | Sweeps **6.89V – 25.39V** |

---

## 3. `C1` Reservoir Capacitor Bench-Test Kit (`470 µF`, `1000 µF`, `1500 µF` @ `25 V`)

All options fit the same **`Capacitor_THT:CP_Radial_D10.0mm_P5.00mm`** (`Ø10 mm`, `5.0 mm` lead pitch) footprint:

| Cap Value | Type | LCSC # | MPN | Brand | Can Size ($\varnothing \times H$) | ESR @ 100kHz | Stock | Price | Notes |
|---|---|---|---|---|---|---:|---:|---:|---|
| **`470 µF 25V`** | Low-Z Electrolytic | **`C108865`** | `25ZLJ470M10X12.5` | Rubycon | `Ø10 × 12.5 mm` | `53 mΩ` | 8,858 | $0.1491 | Compact 12.5 mm height, low-ESR `ZLJ` series |
| **`470 µF 25V`** | Electrolytic | **`C88740`** | `25YXJ470M10X12.5` | Rubycon | `Ø10 × 12.5 mm` | `80 mΩ` | 33,254 | $0.1411 | High stock, compact `Ø10 × 12.5 mm` |
| **`1000 µF 25V`** | **Polymer Aluminum** | **`C46550436`** | `PA25V1000M10x12` | jieerrui | **`Ø10 × 12.5 mm`** | **`16 mΩ`** | 23,232 | $0.4009 | **3.5 mm shorter & lighter** + ultra-low `16 mΩ` ESR |
| **`1000 µF 25V`** | Electrolytic | **`C48935382`** | `ERG25V1000M10X16` | jieerrui | `Ø10 × 16 mm` | `40 mΩ` | 53,755 | $0.1411 | Standard `Ø10 × 16 mm` low-ESR electrolytic |
| **`1000 µF 25V`** | Electrolytic | **`C88720`** | `25YXJ1000M10X20` | Rubycon | `Ø10 × 20 mm` | `46 mΩ` | 48,537 | $0.1995 | Rubycon benchmark (`20 mm` tall) |
| **`1500 µF 25V`** | **Polymer Aluminum** | **`C49232902`** | `PA25V1500M10X15` | jieerrui | **`Ø10 × 15 mm`** | **`16 mΩ`** | 1,834 | $0.5662 | Fits 1500 µF 25V into `15 mm` height with `16 mΩ` ESR |
| **`1500 µF 25V`** | **Polymer Aluminum** | **`C5438719`** | `HT1E158M1016PC` | Honor Elec | **`Ø10 × 16 mm`** | ~`20 mΩ` | 1,451 | $0.4828 | Compact `Ø10 × 16 mm` polymer alternative |
| **`1500 µF 25V`** | Electrolytic | **`C726190`** | `EGF158M1EG20CB` | Man Yue | `Ø10 × 20 mm` | ~`45 mΩ` | 4,040 | $0.1637 | Standard `Ø10 × 20 mm` wet electrolytic |

---

## 4. `R6` Charge-Resistor Bench-Tuning Options (`0805`)

| Value | LCSC # | MPN | Brand | Tier | Stock | Price | Recharge Time (`1000 µF` to 95%) |
|---|---|---|---|---|---:|---:|---|
| **`330 Ω 0805`** *(default)* | **`C17630`** | `0805W8F3300T5E` | UNI-ROYAL | **Basic** | 1,151,059 | $0.0044 | `~1.0 s` (Peak `51 mA` @ 20V — guaranteed below `U1` current limit) |
| **`220 Ω 0805`** *(faster)* | **`C17557`** | `0805W8F2200T5E` | UNI-ROYAL | **Basic** | 747,781 | $0.0029 | `~0.66 s` (Peak `77 mA` @ 20V — near `U1` limit) |
| **`100 Ω 0805`** *(aggressive)* | **`C17408`** | `0805W8F1000T5E` | UNI-ROYAL | **Basic** | 6,284,133 | $0.0034 | Works up to ~14–15V; may trip `U1` 95ms fault >18V on 1000–1500µF |

---

## 5. Refreshing Live Stock

Run from repository root:

```bash
design_optimization/.venv/bin/python scripts/query_jlc_parts.py
```
