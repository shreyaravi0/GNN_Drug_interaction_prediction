Drug Interaction Prediction Project
Overview

This project focuses on predicting drug–drug interactions using a Graph Neural Network (GNN).
Drugs are represented using SMILES (Simplified Molecular Input Line Entry System), which allows the model to analyze molecular structures and learn interaction patterns between drug pairs.

The system combines machine learning, graph-based representations, and deep learning to predict whether two drugs are likely to interact.

drug-interaction-project/
│
├── backend/
│   ├── drug-interaction-aiml.ipynb   # Main notebook for data preprocessing and model training
│   └── main.py                 
|   |__model_def.py
|   |__predictor.py
|   |__utils.py
│
├── dataset/
│   ├── drugbank_train.csv                     # Training dataset containing drug pairs and interaction labels
│   └── drugbank_test.csv                      # Test dataset for model evaluation
│
├── frontend/
└── README.md                         # Project documentation
