import periodictable as pt
from flask import Flask, render_template, request, send_from_directory
from flask_login import login_required, current_user
import os
from element import ELEMENT_DATA
from formula import FORMULAS
from extensions import db, login_manager, bcrypt
from auth import auth

app = Flask(__name__)

# ── Secret key (set this in your Render environment variables too)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'change-this-in-production')

# ── PostgreSQL on Render, SQLite locally
database_url = os.environ.get('DATABASE_URL', '')
if database_url:
    database_url = database_url.replace('postgres://', 'postgresql://')
else:
    database_url = 'sqlite:///chemistrylab.db'  # local fallback

app.config['SQLALCHEMY_DATABASE_URI']        = database_url
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
 
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
def home():
    return render_template("index.html")

@app.route("/periodic-table", methods=["GET", "POST"])
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

@app.route("/virtual")
def virtual_lab():
    return render_template("virtual.html")

@app.route("/contact")
def contact():
    return render_template("contact.html")

@app.route("/exp-acids-bases-salts")
def exp_acids_bases_salts():
    return render_template("exp_abs.html")

@app.route("/exp-lab-safety")
def exp_lab_safety():
    return render_template("exp-lab-safety.html")

@app.route("/exp-titration")
def exp_titration():
    return render_template("exp-titration.html")

if __name__ == "__main__":
    app.run(port=5002, debug=True)
    