import periodictable as pt
from flask import Flask, render_template, request, send_from_directory
from flask_login import login_required, current_user
import os
from element import ELEMENT_DATA
from formula import FORMULAS
from extensions import db, login_manager, bcrypt
from auth import auth
from dotenv import load_dotenv

load_dotenv()

app = Flask(__name__)

# ── Secret key (set this in your Render environment variables too)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'change-this-in-production')

# ── PostgreSQL on Render, SQLite locally
database_url = os.environ.get('DATABASE_URL', '')
if database_url:
    database_url = database_url.replace('postgres://', 'postgresql://')
else:
    database_url = 'sqlite:///chemistrylab.db'  # local fallback


app.config['SQLALCHEMY_DATABASE_URI'] = database_url
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SQLALCHEMY_ENGINE_OPTIONS'] = {
    'pool_pre_ping': True,       # tests connection before using it
    'pool_recycle': 280,         # recycles connections every 280 seconds
    'pool_timeout': 20,          # waits max 20s for a connection
    'pool_size': 5,              # max 5 connections in the pool
    'max_overflow': 2,           # 2 extra connections if pool is full
}
 
# ── Init extensions
db.init_app(app)
bcrypt.init_app(app)
login_manager.init_app(app)
 
# ── Register auth blueprint
app.register_blueprint(auth)
 
# ── Create tables (run once, or use flask db migrate)
with app.app_context():
    db.create_all()
 


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
@login_required
def home():
    return render_template("index.html")

@app.route("/periodic-table", methods=["GET", "POST"])
@login_required
def periodic_table():
    if request.method == "POST":
        name = request.form.get("element", "").strip()
        element = find_element(name)

        if element:
            extra_data = ELEMENT_DATA.get(element.symbol, {})
            return render_template("periodic_table.html", element_data={
                "name": element.name,
                "symbol": element.symbol,
                "number": element.number,
                "density": getattr(element, 'density', None),
                "atomic_mass": element.mass,
                "melting_point": extra_data.get('melting', None),
                "boiling_point": extra_data.get('boiling', None),
                "electronegativity": extra_data.get('electronegativity', None),
                "atomic_radius": extra_data.get('atomic_radius', None),
            })
        else:
            return render_template("invalid.html", source="periodic")

    return render_template("periodic_table.html", element_data=None, source="periodic-table")

@app.route("/notes")
@login_required
def notes():
    return render_template("notes.html")

@app.route("/formula", methods=["GET", "POST"])
@login_required
def formula_search():
    if request.method == "POST":
        query = request.form.get("query", "")
        data = get_formula_data(query)

        if data:
            return render_template("formula.html", data=data)
        else:
            return render_template("invalid.html", source="formula")

    return render_template("formula.html", data=None)

@app.route("/login")
def login():
    return render_template("login.html")

@app.route("/signup")
def signup():
    return render_template("signup.html")

@app.route("/logout")
def logout():
    from flask_login import logout_user
    logout_user()
    return render_template("logout.html")

@app.route("/landing")
def landing():
    return render_template("landing.html")


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
@login_required
def intro_chemistry():
    return render_template("intro.html")

@app.route("/intro-chemistry")
@login_required
def intro():
    return render_template("intro-chemistry.html")

@app.route("/lab-safety")
@login_required
def lab_safety():
    return render_template("lab-safety.html")

@app.route("/scientific-method")
@login_required
def scientific_method():
    return render_template("scientific-method.html")

@app.route("/acids-bases-salts")
@login_required
def acids_bases_and_salts():
    return render_template("abs.html")

@app.route("/lab-equipment")
@login_required
def lab_equipment():
    return render_template("lab-equipment.html")

@app.route("/virtual-lab-guide")
@login_required
def virtual_lab_guide():
    return render_template("virtual-lab-guide.html")

@app.route("/simple-experiments")
@login_required
def simple_experiments():
    return render_template("simple-experiments.html")

@app.route("/observation-recording")
@login_required
def observation_recording():
    return render_template("observation-recording.html")

@app.route("/ss1-chemistry")
@login_required
def ss1_chemistry():
    return render_template("ss1.html")

@app.route("/ss2-chemistry")
@login_required 
def ss2_chemistry():
    return render_template("ss2.html")

@app.route("/ss3-chemistry")
@login_required
def ss3_chemistry():
    return render_template("ss3.html")

@app.route("/chemical-reactions")
@login_required
def chemical_reactions():
    return render_template("chemical-reactions.html")

@app.route("/states-of-matter")
@login_required
def states_of_matter():
    return render_template("states-of-matter.html")

@app.route("/virtual")
@login_required
def virtual_lab():
    return render_template("virtual.html")

@app.route("/contact")
def contact():
    return render_template("contact.html")

@app.route("/exp-acids-bases-salts")
@login_required
def exp_acids_bases_salts():
    return render_template("exp_abs.html")

@app.route("/exp-lab-safety")
@login_required
def exp_lab_safety():
    return render_template("exp-lab-safety.html")

@app.route("/exp-titration")
@login_required
def exp_titrations():
    return render_template("exp_titration.html")

@app.route("/simple-experiments/separating-mixtures")
@login_required
def exp_separating_mixtures():
    return render_template("exp_01.html")

@app.route("/simple-experiments/charles-law")
@login_required
def exp_charles_law():
    return render_template("exp_02.html")

@app.route("/simple-experiments/purity-assessment")
@login_required
def exp_physical_chemical_changes():
    return render_template("exp_03.html")

@app.route("/simple-experiments/simple-distillation")
@login_required 
def exp_simple_distillation():
    return render_template("exp_04.html")

@app.route("/simple-experiments/physical-and-chemical-changes")
@login_required
def exp_double_displacement():
    return render_template("exp_05.html")

@app.route("/simple-experiments/neutralization-reaction")
@login_required
def exp_neutralization_reaction():
    return render_template("exp_06.html")

@app.route("/simple-experiments/displacement-reaction")
@login_required
def exp_heat_of_neutralization():
    return render_template("exp_07.html")

@app.route("/simple-experiments/evaporation")
@login_required
def exp_evaporation():
    return render_template("exp_08.html")

@app.route("/simple-experiments/energy-changes")
@login_required
def exp_titration():
    return render_template("exp_09.html")

@app.route("/simple-experiments/oxygen-gas")
@login_required
def exp_oxygen():
    return render_template("exp_10.html")

@app.route("/simple-experiments/extraction-of-metals")
@login_required
def exp_extraction_of_metal():
    return render_template("exp_11.html")

@app.route("/simple-experiments/electrolysis")
@login_required
def exp_electrolysis():
    return render_template("exp_12.html")

@app.route("/simple-experiments/qualitative-analysis-I")
@login_required
def exp_rate_of_reaction():
    return render_template("exp_13.html")

@app.route("/simple-experiments/qualitative-analysis-II")
@login_required
def exp_qualitative_analysis():
    return render_template("exp_14.html")

@app.route("/simple-experiments/saponification")
@login_required
def exp_saponification():
    return render_template("exp_15.html")

if __name__ == "__main__":
    app.run(port=5002, debug=True)
    