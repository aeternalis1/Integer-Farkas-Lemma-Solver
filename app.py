import os

from flask import Flask
from instance.config import SECRET_KEY

def create_app(test_config=None):
	app = Flask(__name__, instance_relative_config=True)
	app.config.from_pyfile('config.py', silent=True)
	app.secret_key = SECRET_KEY

	import index
	app.register_blueprint(index.bp)

	return app

app = create_app()

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='127.0.0.1', port=port)