default:
	@echo "make run?"

run:
	python3 meshgen_comparison/main.py

clean:
	@find . -name "*.bak" -o -name "__pycache__" -o -name "*.pyc" -o -name "*.pyo" -o -name "*.vtk" -o -name "*.svg" -delete
	@rm -rf *.egg-info/ build/ dist/

format:
	isort .
	black .
	blacken-docs README.md

lint:
	black --check .
	flake8 .
