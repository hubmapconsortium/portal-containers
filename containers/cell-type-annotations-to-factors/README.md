# cell-type-annotations-to-factors

Translate Cell Type Annotation CSV files to Vitessce `factors.json`

Cell Type Annotation CSV file columns:
- cell_id: The cell barcode label from each experiment
- annotation: Our predicted annotation. Each value is a term from the EBI Cell Ontology.
- prediction_score: Confidence level in prediction, ranging from 0..1
