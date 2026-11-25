import json
import urllib.request


# --- Parameters & Search Targets ---

API_URL = "https://jlcpcb.com/api/overseas-pcb-order/v1/shoppingCart/smtGood/selectSmtComponentList"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
    "User-Agent": "Mozilla/5.0",
}

EXISTING_CODES = [
    ("D1", "C2480"),
    ("D2 / D3", "C8598"),
    ("Q1", "C20917"),
    ("R1 (100k 0603)", "C25803"),
    ("R2 (100R 0603)", "C22775"),
    ("C2 (2.2uF 0603)", "C23630"),
    ("R5 (1M Trim Pot)", "C56560"),
    ("L2 (4.7uH Boost)", "C354574"),
    ("Old C3 (22uF 0201)", "C45783"),
    ("Old R3 (30.1k 0201)", "C2907025"),
]

SEARCH_GROUPS = [
    ("U1: TLV61046A Boost Converter (SOT-23-6)", "TLV61046A", None),
    ("R3: 43k 0603 1% Resistor", "43k 0603 1%", "0603"),
    ("R6: 330k 0603 1% Resistor", "330k 0603 1%", "0603"),
    ("R4: 330R 0805 Resistor", "330R 0805 1%", "0805"),
    ("R4 (alt): 220R 0805 Resistor", "220R 0805 1%", "0805"),
    ("R4 (alt): 100R 0805 Resistor", "100R 0805 1%", "0805"),
    ("C3: 10uF 0603 16V/25V X5R/X7R MLCC", "10uF 0603", "0603"),
    ("C3 (alt): 22uF 0603 10V/16V X5R MLCC", "22uF 0603", "0603"),
    ("R5 (1M 3314J Trim Pot): Bourns 3314J 1M", "3314J 1M", None),
    ("R5 (1M 3314J code 105): 3314J-1-105", "3314J-1-105", None),
    ("R4 (if R5=500k): 220k 0603 1% Resistor", "220k 0603 1%", "0603"),
    ("R4 (if R5=500k, 0402): 220k 0402 1% Resistor", "220k 0402 1%", "0402"),
    ("R3 (if R5=500k): 22k 0603 1% Resistor", "22k 0603 1%", "0603"),
    ("R3 (if R5=500k, 0402): 22k 0402 1% Resistor", "22k 0402 1%", "0402"),
    ("R3 (if R5=500k): 24k 0603 1% Resistor", "24k 0603 1%", "0603"),
    ("R3 (if R5=500k, 0402): 24k 0402 1% Resistor", "24k 0402 1%", "0402"),
    ("C4: 4.7uF 0805 50V X5R/X7R MLCC", "4.7uF 0805 50V", "0805"),
    ("C4 (alt): 10uF 0805 35V/50V X5R/X7R MLCC", "10uF 0805 35V", "0805"),
    ("C4 (alt 50V): 10uF 0805 50V X5R MLCC", "10uF 0805 50V", "0805"),
    ("C1 Option A: 470uF 25V Radial Electrolytic", "470uF 25V 10x", None),
    ("C1 Option A2: 470uF 25V Radial (8mm/10mm)", "470uF 25V radial", None),
    ("C1 Option B: 1000uF 25V Radial Electrolytic", "1000uF 25V 10x", None),
    ("C1 Option C: 1500uF 25V Radial Electrolytic", "1500uF 25V", None),
]


def fetch_jlc(keyword, page_size=25):
    payload = json.dumps({
        "currentPage": 1,
        "pageSize": page_size,
        "keyword": keyword,
    }).encode("utf-8")
    req = urllib.request.Request(API_URL, data=payload, headers=HEADERS, method="POST")
    with urllib.request.urlopen(req, timeout=15) as resp:
        body = json.loads(resp.read().decode("utf-8"))
    info = (body.get("data") or {}).get("componentPageInfo") or {}
    return info.get("list") or []


def fmt_item(it):
    code = it.get("componentCode", "")
    mpn = it.get("componentModelEn", "")
    brand = it.get("componentBrandEn", "")
    pkg = it.get("componentSpecificationEn", "")
    lib = it.get("componentLibraryType", "")
    pref = it.get("preferredComponentFlag", False)
    if lib == "base":
        lib_str = "Basic"
    elif pref:
        lib_str = "Preferred"
    else:
        lib_str = "Extended"
    stock = int(it.get("stockCount") or 0)
    prices = it.get("componentPrices") or []
    price_str = f"${prices[0]['productPrice']:.4f}" if prices else "N/A"
    desc = it.get("describe", "")
    return {
        "code": code,
        "mpn": mpn,
        "brand": brand,
        "pkg": pkg,
        "lib": lib_str,
        "stock": stock,
        "price": price_str,
        "desc": desc,
    }


# --- 1. Verify Existing LCSC Codes ---

print("## 1. Existing Schematic LCSC Parts (Live Verification)\n")
print("| Ref | LCSC # | MPN | Brand | Package | JLC Tier | Stock | Unit Price | Key Specs |")
print("|---|---|---|---|---|---|---:|---:|---|")

for ref_label, code in EXISTING_CODES:
    items = fetch_jlc(code, page_size=10)
    exact = [it for it in items if it.get("componentCode") == code]
    if not exact:
        print(f"| {ref_label} | {code} | NOT FOUND | - | - | - | 0 | - | - |")
        continue
    row = fmt_item(exact[0])
    print(
        f"| {ref_label} | {row['code']} | `{row['mpn']}` | {row['brand']} | "
        f"`{row['pkg']}` | {row['lib']} | {row['stock']:,} | {row['price']} | {row['desc']} |"
    )


# --- 2. Search & Rank Candidate Parts ---

print("\n## 2. Candidate Parts Search (In-Stock Ranked by Tier & Stock)\n")

for title, kw, pkg_sub in SEARCH_GROUPS:
    items = fetch_jlc(kw, page_size=30)
    rows = [fmt_item(it) for it in items if int(it.get("stockCount") or 0) > 0]
    if pkg_sub is not None:
        rows = [r for r in rows if pkg_sub in r["pkg"]]

    tier_rank = {"Basic": 0, "Preferred": 1, "Extended": 2}
    rows.sort(key=lambda r: (tier_rank.get(r["lib"], 3), -r["stock"]))

    print(f"### {title} (`{kw}`)")
    print("| LCSC # | MPN | Brand | Package | Tier | Stock | Price | Specs |")
    print("|---|---|---|---|---|---:|---:|---|")
    for r in rows[:5]:
        print(
            f"| {r['code']} | `{r['mpn']}` | {r['brand']} | `{r['pkg']}` | "
            f"{r['lib']} | {r['stock']:,} | {r['price']} | {r['desc']} |"
        )
    print()
