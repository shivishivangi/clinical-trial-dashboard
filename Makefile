setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis.py

dashboard:
	python dashboard.py

test:
	python -m unittest discover -s tests -p "test_*.py" -v

.PHONY: setup pipeline dashboard test