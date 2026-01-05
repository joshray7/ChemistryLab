import periodictable as pt
from flask import Flask, render_template, request, send_from_directory
import os

app = Flask(__name__)

FORMULAS = {
    "ohms law": {
        "formula": "V = I × R",
        "symbols": {
            "V": "Voltage (volts)",
            "I": "Current (amperes)",
            "R": "Resistance (ohms)"
        }
    },

    "boyle's law": {
        "formula": "P₁V₁ = P₂V₂",
        "symbols": {
            "P": "Pressure",
            "V": "Volume",
            "1": "Initial state",
            "2": "Final state"
        }
    },

    "charles law": {
        "formula": "V₁ / T₁ = V₂ / T₂",
        "symbols": {
            "V": "Volume",
            "T": "Temperature (Kelvin)",
            "1": "Initial state",
            "2": "Final state"
        }
    },

    "avogadro's law": {
        "formula": "V₁ / n₁ = V₂ / n₂",
        "symbols": {
            "V": "Volume",
            "n": "Number of moles",
            "1": "Initial state",
            "2": "Final state"
        }
    },

    "gay lussac law": {
        "formula": "P₁ / T₁ = P₂ / T₂",
        "symbols": {
            "P": "Pressure",
            "T": "Temperature (Kelvin)",
            "1": "Initial state",
            "2": "Final state"
        }
    },

    "combined gas law": {
        "formula": "P₁V₁ / T₁ = P₂V₂ / T₂",
        "symbols": {
            "P": "Pressure",
            "V": "Volume",
            "T": "Temperature"
        }
    },

    "ideal gas law": {
        "formula": "PV = nRT",
        "symbols": {
            "P": "Pressure",
            "V": "Volume",
            "n": "Number of moles",
            "R": "Gas constant",
            "T": "Temperature (Kelvin)"
        }
    },

    "dalton's law": {
        "formula": "P_total = P₁ + P₂ + ...",
        "symbols": {
            "P_total": "Total pressure",
            "P₁": "Partial pressure of gas 1",
            "P₂": "Partial pressure of gas 2"
        }
    },

    "graham's law": {
        "formula": "r₁ / r₂ = √(M₂ / M₁)",
        "symbols": {
            "r": "Rate of diffusion",
            "M": "Molar mass"
        }
    },

    "henry's law": {
        "formula": "C = kP",
        "symbols": {
            "C": "Concentration",
            "k": "Henry's constant",
            "P": "Pressure"
        }
    },

    "molarity": {
        "formula": "M = mol / L",
        "symbols": {
            "M": "Molarity",
            "mol": "Number of moles",
            "L": "Volume (liters)"
        }
    },

    "molality": {
        "formula": "m = mol / kg",
        "symbols": {
            "m": "Molality",
            "mol": "Moles of solute",
            "kg": "Mass of solvent (kg)"
        }
    },

    "normality": {
        "formula": "N = equivalents / L",
        "symbols": {
            "N": "Normality",
            "equivalents": "Equivalent mass",
            "L": "Volume (liters)"
        }
    },

    "dilution formula": {
        "formula": "M₁V₁ = M₂V₂",
        "symbols": {
            "M": "Molarity",
            "V": "Volume",
            "1": "Initial state",
            "2": "Final state"
        }
    },

    "percentage composition": {
        "formula": "% = (mass of solute / mass of solution) × 100",
        "symbols": {
            "mass of solute": "Mass of solute",
            "mass of solution": "Total mass"
        }
    },

    "ppm": {
        "formula": "ppm = (mass solute / mass solution) × 10⁶",
        "symbols": {}
    },

    "ppb": {
        "formula": "ppb = (mass solute / mass solution) × 10⁹",
        "symbols": {}
    },

    "osmotic pressure": {
        "formula": "π = MRT",
        "symbols": {
            "π": "Osmotic pressure",
            "M": "Molarity",
            "R": "Gas constant",
            "T": "Temperature"
        }
    },

    "freezing point depression": {
        "formula": "ΔT_f = K_f m",
        "symbols": {
            "ΔT_f": "Change in freezing point",
            "K_f": "Cryoscopic constant",
            "m": "Molality"
        }
    },

    "heat equation": {
        "formula": "Q = mcΔT",
        "symbols": {
            "Q": "Heat energy",
            "m": "Mass",
            "c": "Specific heat capacity",
            "ΔT": "Temperature change"
        }
    },

    "gibb's free energy": {
        "formula": "ΔG = ΔH − TΔS",
        "symbols": {
            "ΔG": "Gibbs free energy",
            "ΔH": "Enthalpy change",
            "T": "Temperature",
            "ΔS": "Entropy change"
        }
    },

    "work done": {
        "formula": "W = PΔV",
        "symbols": {
            "W": "Work",
            "P": "Pressure",
            "ΔV": "Change in volume"
        }
    },

    "equilibrium constant": {
        "formula": "Kc = [products] / [reactants]",
        "symbols": {
            "Kc": "Equilibrium constant"
        }
    },

    "rate law": {
        "formula": "Rate = k[A]^m[B]^n",
        "symbols": {
            "k": "Rate constant",
            "A": "Concentration of A",
            "B": "Concentration of B",
            "m": "Order for A",
            "n": "Order for B"
        }
    },

    "arrhenius equation": {
        "formula": "k = Ae^(−Ea/RT)",
        "symbols": {
            "k": "Rate constant",
            "A": "Frequency factor",
            "Ea": "Activation energy",
            "R": "Gas constant",
            "T": "Temperature"
        }
    },

    "combined gas law": {
        "formula": "P₁V₁ / T₁ = P₂V₂ / T₂",
        "symbols": {
            "P": "Pressure",
            "V": "Volume",
            "T": "Temperature (Kelvin)",
            "1": "Initial state",
            "2": "Final state"
        }
    },

    "ideal gas law": {
        "formula": "PV = nRT",
        "symbols": {
            "P": "Pressure",
            "V": "Volume",
            "n": "Number of moles",
            "R": "Gas constant",
            "T": "Temperature (Kelvin)"
        }
    },

    "newton's first law": {
        "formula": "F = 0 (when ΣF = 0)",
        "symbols": {
            "F": "Force",
            "ΣF": "Resultant force"
        }
    },

    "newton's second law": {
        "formula": "F = m × a",
        "symbols": {
            "F": "Force (newtons)",
            "m": "Mass (kg)",
            "a": "Acceleration (m/s²)"
        }
    },

    "newton's third law": {
        "formula": "F₁ = -F₂",
        "symbols": {
            "F₁": "Action force",
            "F₂": "Reaction force"
        }
    },

    "weight": {
        "formula": "W = m × g",
        "symbols": {
            "W": "Weight (newtons)",
            "m": "Mass (kg)",
            "g": "Acceleration due to gravity (m/s²)"
        }
    },

    "density": {
        "formula": "ρ = m / V",
        "symbols": {
            "ρ": "Density (kg/m³)",
            "m": "Mass (kg)",
            "V": "Volume (m³)"
        }
    },

    "pressure": {
        "formula": "P = F / A",
        "symbols": {
            "P": "Pressure (N/m²)",
            "F": "Force (newtons)",
            "A": "Area (m²)"
        }
    },

    "moment of a force": {
        "formula": "M = F × d",
        "symbols": {
            "M": "Moment (Nm)",
            "F": "Force (newtons)",
            "d": "Perpendicular distance (m)"
        }
    },

    "kinetic energy": {
        "formula": "KE = 1/2 m v²",
        "symbols": {
            "KE": "Kinetic energy (joules)",
            "m": "Mass (kg)",
            "v": "Velocity (m/s)"
        }
    },

    "potential energy": {
        "formula": "PE = m g h",
        "symbols": {
            "PE": "Potential energy (joules)",
            "m": "Mass (kg)",
            "g": "Acceleration due to gravity",
            "h": "Height (m)"
        }
    },

    "work done": {
        "formula": "W = F × d",
        "symbols": {
            "W": "Work done (joules)",
            "F": "Force (newtons)",
            "d": "Distance moved (m)"
        }
    },

    "power": {
        "formula": "P = W / t",
        "symbols": {
            "P": "Power (watts)",
            "W": "Work done (joules)",
            "t": "Time (seconds)"
        }
    },

    "speed": {
        "formula": "v = d / t",
        "symbols": {
            "v": "Speed (m/s)",
            "d": "Distance (m)",
            "t": "Time (s)"
        }
    },

    "velocity": {
        "formula": "v = Δs / Δt",
        "symbols": {
            "v": "Velocity (m/s)",
            "Δs": "Displacement",
            "Δt": "Time interval"
        }
    },

    "acceleration": {
        "formula": "a = Δv / Δt",
        "symbols": {
            "a": "Acceleration (m/s²)",
            "Δv": "Change in velocity",
            "Δt": "Time interval"
        }
    },

    "ohmic power law": {
        "formula": "P = V × I",
        "symbols": {
            "P": "Power (watts)",
            "V": "Voltage (volts)",
            "I": "Current (amperes)"
        }
    },

    "energy consumption": {
        "formula": "E = P × t",
        "symbols": {
            "E": "Energy (joules)",
            "P": "Power (watts)",
            "t": "Time (seconds)"
        }
    },

    "snell's law": {
        "formula": "n₁ sinθ₁ = n₂ sinθ₂",
        "symbols": {
            "n₁": "Refractive index of medium 1",
            "n₂": "Refractive index of medium 2",
            "θ₁": "Angle of incidence",
            "θ₂": "Angle of refraction"
        }
    },

    "magnification lens": {
        "formula": "m = hᵢ / hₒ = v / u",
        "symbols": {
            "m": "Magnification",
            "hᵢ": "Image height",
            "hₒ": "Object height",
            "v": "Image distance",
            "u": "Object distance"
        }
    },

    "lens formula": {
        "formula": "1/f = 1/v + 1/u",
        "symbols": {
            "f": "Focal length",
            "v": "Image distance",
            "u": "Object distance"
        }
    },

    "mirror formula": {
        "formula": "1/f = 1/v + 1/u",
        "symbols": {
            "f": "Focal length",
            "v": "Image distance",
            "u": "Object distance"
        }
    },

    "refractive index": {
        "formula": "n = c / v",
        "symbols": {
            "n": "Refractive index",
            "c": "Speed of light in vacuum",
            "v": "Speed of light in medium"
        }
    },

    "hooke's law": {
        "formula": "F = kx",
        "symbols": {
            "F": "Force (newtons)",
            "k": "Spring constant",
            "x": "Extension or compression (m)"
        }
    },

    "stress": {
        "formula": "σ = F / A",
        "symbols": {
            "σ": "Stress (N/m²)",
            "F": "Force",
            "A": "Cross-sectional area"
        }
    },

    "strain": {
        "formula": "ε = ΔL / L",
        "symbols": {
            "ε": "Strain",
            "ΔL": "Change in length",
            "L": "Original length"
        }
    },

    "young's modulus": {
        "formula": "E = σ / ε",
        "symbols": {
            "E": "Young's modulus",
            "σ": "Stress",
            "ε": "Strain"
        }
    },

    "coulombs law": {
        "formula": "F = k(q₁q₂) / r²",
        "symbols": {
            "F": "Electric force (newtons)",
            "k": "Coulomb's constant",
            "q₁": "Charge 1",
            "q₂": "Charge 2",
            "r": "Distance between charges"
        }
    },

    "electrical resistance": {
        "formula": "R = ρL / A",
        "symbols": {
            "R": "Resistance (ohms)",
            "ρ": "Resistivity",
            "L": "Length of conductor",
            "A": "Cross-sectional area"
        }
    },

    "capacitance": {
        "formula": "C = Q / V",
        "symbols": {
            "C": "Capacitance (farads)",
            "Q": "Charge (coulombs)",
            "V": "Voltage (volts)"
        }
    },

    "charging capacitor voltage": {
        "formula": "V = V₀(1 - e^{-t/RC})",
        "symbols": {
            "V": "Voltage at time t",
            "V₀": "Maximum voltage",
            "R": "Resistance",
            "C": "Capacitance",
            "t": "Time"
        }
    },

    "frequency": {
        "formula": "f = 1 / T",
        "symbols": {
            "f": "Frequency (Hz)",
            "T": "Time period (s)"
        }
    },

    "wave speed": {
        "formula": "v = fλ",
        "symbols": {
            "v": "Wave speed",
            "f": "Frequency",
            "λ": "Wavelength"
        }
    },

    "momentum": {
        "formula": "p = m × v",
        "symbols": {
            "p": "Momentum",
            "m": "Mass",
            "v": "Velocity"
        }
    },

    "impulse": {
        "formula": "I = F × t",
        "symbols": {
            "I": "Impulse",
            "F": "Force",
            "t": "Time"
        }
    },

    "heat energy": {
        "formula": "Q = mcΔT",
        "symbols": {
            "Q": "Heat energy",
            "m": "Mass",
            "c": "Specific heat capacity",
            "ΔT": "Change in temperature"
        }
    },

    "latent heat": {
        "formula": "Q = mL",
        "symbols": {
            "Q": "Heat energy",
            "m": "Mass",
            "L": "Specific latent heat"
        }
    },

    "efficiency": {
        "formula": "Efficiency = (Useful energy output / Total energy input) × 100%",
        "symbols": {
            "Useful energy output": "Energy used for useful work",
            "Total energy input": "Energy supplied"
        }
    },

    "density": {
        "formula": "ρ = m / V",
        "symbols": {
            "ρ": "Density",
            "m": "Mass",
            "V": "Volume"
        }
    },

    "pressure": {
        "formula": "P = F / A",
        "symbols": {
            "P": "Pressure",
            "F": "Force",
            "A": "Area"
        }
    },

    "work": {
        "formula": "W = F × d",
        "symbols": {
            "W": "Work done (joules)",
            "F": "Force",
            "d": "Distance moved"
        }
    },

    "power": {
        "formula": "P = W / t",
        "symbols": {
            "P": "Power (watts)",
            "W": "Work done",
            "t": "Time"
        }
    },

    "kinetic energy": {
        "formula": "KE = 1/2 m v²",
        "symbols": {
            "KE": "Kinetic energy",
            "m": "Mass",
            "v": "Velocity"
        }
    },

    "potential energy": {
        "formula": "PE = m g h",
        "symbols": {
            "PE": "Potential energy",
            "m": "Mass",
            "g": "Acceleration due to gravity",
            "h": "Height"
        }
    },

    "specific heat capacity": {
        "formula": "c = Q / (m ΔT)",
        "symbols": {
            "c": "Specific heat capacity",
            "Q": "Heat energy",
            "m": "Mass",
            "ΔT": "Change in temperature"
        }
    },

    "gas constant": {
        "formula": "PV = nRT",
        "symbols": {
            "P": "Pressure",
            "V": "Volume",
            "n": "Moles of gas",
            "R": "Gas constant",
            "T": "Temperature"
        }
    },

    "ph calculation": {
        "formula": "pH = -log[H⁺]",
        "symbols": {
            "[H⁺]": "Hydrogen ion concentration"
        }
    },

    "poh calculation": {
        "formula": "pOH = -log[OH⁻]",
        "symbols": {
            "[OH⁻]": "Hydroxide ion concentration"
        }
    },

    "ionic product of water": {
        "formula": "Kw = [H⁺][OH⁻]",
        "symbols": {
            "Kw": "Ionic product of water",
            "[H⁺]": "Hydrogen ion concentration",
            "[OH⁻]": "Hydroxide ion concentration"
        }
    },

    "ideal gas law": {
        "formula": "PV = nRT",
        "symbols": {
            "P": "Pressure",
            "V": "Volume",
            "n": "Moles",
            "R": "Gas constant",
            "T": "Temperature (Kelvin)"
        }
    },

    "dilution formula": {
        "formula": "M₁V₁ = M₂V₂",
        "symbols": {
            "M₁": "Initial molarity",
            "V₁": "Initial volume",
            "M₂": "Final molarity",
            "V₂": "Final volume"
        }
    },

    "molarity": {
        "formula": "M = n / V",
        "symbols": {
            "M": "Molarity",
            "n": "Moles of solute",
            "V": "Volume of solution (L)"
        }
    },

    "molality": {
        "formula": "m = n / kg solvent",
        "symbols": {
            "m": "Molality",
            "n": "Moles of solute",
            "kg solvent": "Kilograms of solvent"
        }
    },

    "percent yield": {
        "formula": "Percent yield = (Actual yield / Theoretical yield) × 100",
        "symbols": {
            "Actual yield": "Amount produced",
            "Theoretical yield": "Maximum possible amount"
        }
    },

    "boyle's constant form": {
        "formula": "P ∝ 1/V",
        "symbols": {
            "P": "Pressure",
            "V": "Volume"
        }
    },

    "charles constant form": {
        "formula": "V ∝ T",
        "symbols": {
            "V": "Volume",
            "T": "Temperature"
        }
    },

    "avogadro's constant form": {
        "formula": "V ∝ n",
        "symbols": {
            "V": "Volume",
            "n": "Moles of gas"
        }
    },

    "faraday's first law": {
        "formula": "m = (Q × M) / (F × z)",
        "symbols": {
            "m": "Mass deposited",
            "Q": "Total charge",
            "M": "Molar mass",
            "F": "Faraday constant",
            "z": "Number of electrons transferred"
        }
    },


    "faradays second law": {
        "formula": "m₁ / m₂ = M₁ / M₂",
        "symbols": {
            "m₁": "Mass of first substance",
            "m₂": "Mass of second substance",
            "M₁": "Molar mass of first substance",
            "M₂": "Molar mass of second substance"
        }
    },

    "coulombs law": {
        "formula": "F = k(q₁q₂ / r²)",
        "symbols": {
            "F": "Force between charges",
            "k": "Coulomb's constant",
            "q₁": "Charge 1",
            "q₂": "Charge 2",
            "r": "Distance between charges"
        }
    },

    "frequency": {
        "formula": "f = 1 / T",
        "symbols": {
            "f": "Frequency",
            "T": "Time period"
        }
    },

    "speed": {
        "formula": "v = d / t",
        "symbols": {
            "v": "Speed",
            "d": "Distance",
            "t": "Time"
        }
    },

    "newton's second law": {
        "formula": "F = m a",
        "symbols": {
            "F": "Force",
            "m": "Mass",
            "a": "Acceleration"
        }
    },

    "moment": {
        "formula": "M = F × d",
        "symbols": {
            "M": "Moment",
            "F": "Force",
            "d": "Perpendicular distance"
        }
    },

    "pressure in liquids": {
        "formula": "P = ρ g h",
        "symbols": {
            "P": "Pressure",
            "ρ": "Density",
            "g": "Acceleration due to gravity",
            "h": "Height/depth"
        }
    },

    "ohmic power": {
        "formula": "P = I² R",
        "symbols": {
            "P": "Power",
            "I": "Current",
            "R": "Resistance"
        }
    },

    "electrical energy": {
        "formula": "E = P t",
        "symbols": {
            "E": "Electrical energy",
            "P": "Power",
            "t": "Time"
        }
    },

    "gravitational force": {
        "formula": "F = (G m₁ m₂) / r²",
        "symbols": {
            "F": "Gravitational force",
            "G": "Gravitational constant",
            "m₁": "Mass 1",
            "m₂": "Mass 2",
            "r": "Distance between masses"
        }
    },

    "wave speed": {
        "formula": "v = f λ",
        "symbols": {
            "v": "Wave speed",
            "f": "Frequency",
            "λ": "Wavelength"
        }
    },

    "refractive index": {
        "formula": "n = sin i / sin r",
        "symbols": {
            "n": "Refractive index",
            "i": "Angle of incidence",
            "r": "Angle of refraction"
        }
    },

    "efficiency mechanical": {
        "formula": "Efficiency = (Useful energy output / Total energy input) × 100",
        "symbols": {
            "Useful energy output": "Energy used for work",
            "Total energy input": "All energy supplied"
        }
    },

    "efficiency machines": {
        "formula": "Efficiency = (MA / VR) × 100",
        "symbols": {
            "MA": "Mechanical advantage",
            "VR": "Velocity ratio"
        }
    },

    "mechanical advantage": {
        "formula": "MA = Load / Effort",
        "symbols": {
            "Load": "Weight lifted",
            "Effort": "Force applied"
        }
    },

    "velocity ratio": {
        "formula": "VR = Distance moved by effort / Distance moved by load",
        "symbols": {
            "Distance moved by effort": "Effort movement",
            "Distance moved by load": "Load movement"
        }
    },

    "hooke's law": {
        "formula": "F = k x",
        "symbols": {
            "F": "Force applied",
            "k": "Spring constant",
            "x": "Extension"
        }
    },

    "heat capacity": {
        "formula": "C = Q / ΔT",
        "symbols": {
            "C": "Heat capacity",
            "Q": "Heat energy",
            "ΔT": "Temperature change"
        }
    },

    "charge": {
        "formula": "Q = I t",
        "symbols": {
            "Q": "Charge",
            "I": "Current",
            "t": "Time"
        }
    },

    "universal gas equation per mole": {
        "formula": "PV = NkT",
        "symbols": {
            "P": "Pressure",
            "V": "Volume",
            "N": "Number of molecules",
            "k": "Boltzmann constant",
            "T": "Temperature"
        }
    }

}

def get_formula_data(query):
    query = query.lower().strip()

    # Direct match (best)
    if query in FORMULAS:
        return FORMULAS[query]

    # Try partial match (e.g. user types "boyle" instead of "boyles law")
    for key in FORMULAS:
        if query in key:
            return FORMULAS[key]

    return None


def find_element(name):
    name = name.strip()

    # Try search by symbol (H, O, Na...)
    for el in pt.elements:
        if el.symbol.lower() == name.lower():
            return el

    # Try search by element name (Hydrogen, Oxygen...)
    for el in pt.elements:
        if el.name.lower() == name.lower():
            return el

    return None   # If not found


def clean_formula_response(raw):
    # Split around | (Wolfram uses | as separators)
    parts = raw.split("|")
    
    # First part = the formula
    formula = parts[0].strip()
    
    definitions = []

    for part in raw.split("\n"):
        # If Wolfram gives "V | voltage", "I | current", etc.
        if "|" in part and len(part.split("|")) == 2:
            definitions.append(part.strip())


    return formula, definitions


@app.route("/google6bd543ba6834f0bf.html")
def google_verification():
    return send_from_directory(os.getcwd(), 'google6bd543ba6834f0bf.html')

@app.route("/")
def home():
    return render_template("index.html")


@app.route("/periodic-table", methods=["GET", "POST"])
def periodic_table():
    if request.method == "POST":
        name = request.form.get("element", "").strip()
        element = find_element(name)

        if element:
            return render_template("periodic_table.html", element_data={
                "name": element.name,
                "symbol": element.symbol,
                "number": element.number,
                "density": element.density,
            })
        else:
            return render_template("invalid.html", source="periodic")

    return render_template("periodic_table.html", element_data=None, source="periodic-table")

@app.route("/notes")
def notes():
    return render_template("notes.html")

@app.route("/formula", methods=["GET", "POST"])
def formula_search():
    if request.method == "POST":
        query = request.form.get("query", "")
        data = get_formula_data(query)

        if data:
            return render_template("formula.html", data=data)
        else:
            return render_template("invalid.html", source="formula")

    return render_template("formula.html", data=None)

@app.route("/feedback", methods=["GET", "POST"])
def feedback():
    return render_template("feedback.html")

@app.route("/about-us", methods=["GET", "POST"])
def about():
    return render_template("about.html")

@app.route("/privacy-policy")
def privacy_policy():
    return render_template("privacy-policy.html")

@app.route("/terms")
def terms():
    return render_template("terms.html")

@app.route("/intro")
def intro_chemistry():
    return render_template("intro.html")

@app.route("/intro-chemistry")
def intro():
    return render_template("intro-chemistry.html")

@app.route("/lab-safety")
def lab_safety():
    return render_template("lab-safety.html")

@app.route("/scientific-method")
def scientific_method():
    return render_template("scientific-method.html")

@app.route("/acids-bases-salts")
def acids_bases_and_salts():
    return render_template("abs.html")

@app.route("/lab-equipment")
def lab_equipment():
    return render_template("lab-equipment.html")

@app.route("/virtual-lab-guide")
def virtual_lab_guide():
    return render_template("virtual-lab-guide.html")

@app.route("/simple-experiments")
def simple_experiments():
    return render_template("simple-experiments.html")

@app.route("/observation-recording")
def observation_recording():
    return render_template("observation-recording.html")

@app.route("/ss1-chemistry")
def ss1_chemistry():
    return render_template("ss1.html")

@app.route("/ss2-chemistry")
def ss2_chemistry():
    return render_template("ss2.html")

@app.route("/ss3-chemistry")
def ss3_chemistry():
    return render_template("ss3.html")

@app.route("/chemical-reactions")
def chemical_reactions():
    return render_template("chemical-reactions.html")

@app.route("/states-of-matter")
def states_of_matter():
    return render_template("states-of-matter.html")

if __name__ == "__main__":
    app.run(debug=True)
    