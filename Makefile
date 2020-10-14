clean:
	@find . -name "*.bak" -o -name "__pycache__" -o -name "*.pyc" -o -name "*.pyo" -o -name "*.vtk" -delete
	@rm -rf *.egg-info/ build/ dist/

format:
	isort .
	black .
	blacken-docs README.md

lint:
	black --check .
	flake8 .
