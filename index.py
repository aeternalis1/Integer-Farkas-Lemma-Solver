from flask import (
	Blueprint, flash, g, redirect, render_template, request, session, url_for
)

from flask_wtf import FlaskForm
from wtforms import IntegerField, StringField, SubmitField
from wtforms.validators import DataRequired, Length
from wtforms.widgets import TextArea

from scripts.parse import ParseInput, InputError
from scripts.main import Solve

bp = Blueprint('index', __name__)


class MatrixForm(FlaskForm):
    """Contact form."""
    m = IntegerField(
        'Number of rows (m)',
        [DataRequired()]
    )
    n = IntegerField(
        'Number of columns (n)',
        [DataRequired()]
    )
    mat = StringField(
        'A rational matrix of dimensions m x (n+1)',
        [DataRequired()],
        widget=TextArea()
    )
    submit = SubmitField('Submit')


@bp.route('/', methods=['GET', 'POST'])
def index():
	form = MatrixForm()

	if form.validate_on_submit():
		m = form.m.data
		n = form.n.data
		mat = form.mat.data 
		try:
			M = ParseInput(m, n, mat)
			res = Solve(M)
			print ("res:", res.t)
			return render_template("result.html", m=form.m.data, n=form.n.data, res=res)
		except InputError as e:
			flash("Invalid input: %s" % e)
	return render_template(
		"index.html",
		form=form,
		template="form-template"
	)
