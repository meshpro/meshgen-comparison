clean:
	@find . | grep -E "(__pycache__|\.pyc|\.pyo|\.vtu$\)" | xargs rm -rf
	@rm -rf *.egg-info/ build/ dist/

format:
	isort .
	black .
	blacken-docs README.md

lint:
	black --check .
	flake8 .
