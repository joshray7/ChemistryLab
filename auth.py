from flask import Blueprint, render_template, redirect, url_for, flash, request
from flask_login import login_user, logout_user, login_required, current_user
from extensions import db, bcrypt
from models import User
 
auth = Blueprint('auth', __name__)
 
 
# ── SIGNUP ──────────────────────────────────────────────────
@auth.route('/signup', methods=['GET', 'POST'])
def signup():
    if current_user.is_authenticated:
        return redirect(url_for('home'))
 
    if request.method == 'POST':
        username = request.form.get('username', '').strip()
        password = request.form.get('password', '').strip()
        confirm  = request.form.get('confirm_password', '').strip()
 
        # ── Validation
        error = None
        if not username or not password or not confirm:
            error = 'All fields are required.'
        elif len(username) < 3:
            error = 'Username must be at least 3 characters.'
        elif len(username) > 80:
            error = 'Username must be under 80 characters.'
        elif len(password) < 6:
            error = 'Password must be at least 6 characters.'
        elif password != confirm:
            error = 'Passwords do not match.'
        elif User.query.filter_by(username=username).first():
            error = 'That username is already taken.'
 
        if error:
            flash(error, 'error')
            return render_template('signup.html', username=username)
 
        # ── Create user
        hashed_pw = bcrypt.generate_password_hash(password).decode('utf-8')
        new_user  = User(username=username, password=hashed_pw)
        db.session.add(new_user)
        db.session.commit()
 
        flash('Account created! You can now log in.', 'success')
        return redirect(url_for('auth.login'))
 
    return render_template('signup.html')
 
 
# ── LOGIN ────────────────────────────────────────────────────
@auth.route('/login', methods=['GET', 'POST'])
def login():
    if current_user.is_authenticated:
        return redirect(url_for('home'))
 
    if request.method == 'POST':
        username = request.form.get('username', '').strip()
        password = request.form.get('password', '').strip()
        remember = request.form.get('remember') == 'on'
 
        error = None
        if not username or not password:
            error = 'Both fields are required.'
 
        if not error:
            user = User.query.filter_by(username=username).first()
            if not user or not bcrypt.check_password_hash(user.password, password):
                error = 'Incorrect username or password.'
 
        if error:
            flash(error, 'error')
            return render_template('login.html', username=username)
 
        login_user(user, remember=remember)
        next_page = request.args.get('next')
        if current_user.is_authenticated:
            return redirect(url_for('home'))
 
    return render_template('login.html')
 
 
# ── LOGOUT ───────────────────────────────────────────────────
@auth.route('/logout')
@login_required
def logout():
    logout_user()
    flash('You have been logged out.', 'info')
    return redirect(url_for('auth.login'))
