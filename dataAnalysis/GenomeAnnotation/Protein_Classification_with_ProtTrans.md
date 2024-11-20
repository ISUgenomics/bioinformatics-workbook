---
title: A ProtTrans Pipeline to Differentiate Transmembrane Domains from Other Proteins
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---



# A ProtTrans pipeline to differentiate transmembrane domains containing proteins from others
### Introduction to ProtTrans for Bioinformatics Applications

ProtTrans is a powerful tool that combines cutting-edge artificial intelligence (AI) with protein sequence analysis. ProtTrans enables bioinformaticians to analyze and predict protein functions, structures, and classifications in ways that were previously unattainable. Whether you're working with kinases, functional domains, or even evolutionary questions, ProtTrans opens up new possibilities.

### What is ProtTrans?

ProtTrans is a transformer-based deep learning model designed to analyze protein sequences. It’s inspired by natural language processing (NLP) tools like BERT (Bidirectional Encoder Representations from Transformers) and adapts them for biological sequences. The model is trained on massive datasets of protein sequences from databases like UniProt, learning patterns, relationships, and evolutionary signals encoded in the amino acid sequences. These learned patterns are used to generate embeddings—numerical representations of proteins that capture their functional and structural properties.

### How ProtTrans Works and Where AI Fits in

**Deep Learning for Embeddings**
*  ProtTrans uses transformers, a type of deep learning AI, to generate protein embeddings.
    *  This step is analogous to converting raw text into dense numerical representations in NLP. Here, proteins are the "text," and the embeddings capture biochemical and functional properties.
*  Machine Learning for Downstream Tasks:
    *  After generating embeddings, you can feed them into traditional machine learning algorithms like logistic regression, random forests, or neural networks for classification, clustering, or prediction tasks.
    *  Example: Classifying proteins as transmembrane or non-transmembrane using a supervised model.
*  Data Preprocessing:
    *  AI-based preprocessing, like scaling and dimensionality reduction, ensures the embeddings are suitable for downstream machine learning models.

Each step combines insights from AI with domain knowledge, enabling you to tackle complex biological questions effectively.

### Applications of ProtTrans

*  ProtTrans has a wide range of applications in bioinformatics, including but not limited to:

    *  Protein Function Prediction:
        *  Determine the function of uncharacterized proteins by analyzing their sequence embeddings.

    *  Protein Classification:
        *  Classify proteins into functional categories, such as kinases, receptors, or transporters.

    *  Protein Structure and Folding:
        *  Leverage embeddings to predict structural features like secondary structure or intrinsic disorder.

    *  Variant Effect Prediction:
        *  Study the effects of mutations by comparing embeddings of wild-type and mutant sequences.

    *  Drug Discovery and Target Identification:
        *  Identify potential therapeutic targets by analyzing protein families, domains, or functional sites.

    *  Comparative Genomics:
        *  Use embeddings for clustering and similarity analysis across species to study protein evolution and conservation.

    *  Pathway Analysis:
        *  Annotate proteins involved in specific biochemical pathways by integrating ProtTrans outputs with pathway databases.

### Why Use ProtTrans

*  Scalability:
    *  Analyze thousands to millions of sequences efficiently.

*  Flexibility:
    *  Apply embeddings to diverse machine learning workflows for predictions, clustering, or annotations.

*  State-of-the-Art AI:
    *  ProtTrans leverages deep learning, ensuring high accuracy and biological relevance.

*  Biological Insights:
    *  Embeddings capture complex relationships in protein sequences that go beyond traditional sequence alignment methods.


### Learning Goals

*  Understand how ProtTrans uses AI to generate protein embeddings.
*  Learn how to integrate these embeddings into machine learning workflows for specific biological questions.
*  Gain insights into the versatility of ProtTrans for bioinformatics research.


### Installation of ProtTrans
```
#create a python virtual environment
python3 -m venv ProtTrans_pyenv
#activate the environment
source ProtTrans_pyenv/bin/activate


#install these three packages
pip install --upgrade pip
pip install -q transformers
pip install torch
pip install sentencepiece
pip install protobuf
pip install h5py
pip install pandas
pip install sklearn
pip install scikit-learn
```

### Create a dataset for training

Here I just went to Uniprot.org to find proteins that were representative transmembrane domain containing proteins from UniRef (manually reviewed). Since we have been working with Arabidopsis, I have kept with that lineage by extracting only proteins from the Brassicaceae. To avoid biasing the model, we want to keep both protein types at an equal representation. 


### Transmembrane dataset
```

# number of proteins in each file. 
grep -c ">" *ransmembrane*.fasta
NotTransmembrane3700.fasta:64601
Transmembrane3700.fasta:18383

#The transmembrane dataset had proteins without transmembrane domains.  The nontransmembrane dataset had proteins with transmembrane domains. 

awk '{print $1}' ../NotTransmembrane3700.fasta >CleanTransmembrane3700.fasta
awk '{print $1}' ../NotTransmembrane3700.fasta >CleanNonTransmembrane3700.fasta
echo "tmhmm CleanNonTransmembrane3700.fasta" >NotTM.sh
echo "tmhmm CleanTransmembrane3700.fasta" >>NotTM.sh

grep "Number of predicted TMHs" CleanNonTransmembrane3700.tmhmmout |awk '$7==0 {print $2}' |cdbyank CleanNonTransmembrane3700.fasta.cidx |awk '{print $1}' >ActualTransmembraneNegative.fasta

grep "Number of predicted TMHs" CleanTransmembrane3700.tmhmmout |awk '$7!=0 {print $2}' |cdbyank CleanTransmembrane3700.fasta.cidx |awk '{print $1}' >ActualTransmembranePositive.fasta


# the distribution of proteins from each class for training
grep -c ">" *fasta
ActualTransmembraneNegative.fasta:14997
ActualTransmembranePositive.fasta:15352

cat ActualTransmembranePositive.fasta ActualTransmembraneNegative.fasta >TMProteinDataset.fasta 
```


# Embedding Protein Sequences

### Prepare sequences for embedding the best protein model provided
**Copy this script to a file named: generate_embeddings.py**
```
import torch
from transformers import AutoTokenizer, AutoModel
import pickle
import argparse

# Input Parameters
parser = argparse.ArgumentParser(description="Generate ProtTrans embeddings")
parser.add_argument('--input', type=str, required=True, help="Input FASTA file")
parser.add_argument('--output', type=str, required=True, help="Output pickle file for embeddings")
args = parser.parse_args()


#Here the tokenizer is converting the protein sequence into tokens, with each token representing a single amino acid using prot_bert_bfd, a deep learning transformer model.
def generate_embeddings(input_fasta, output_pkl):
    tokenizer = AutoTokenizer.from_pretrained("Rostlab/prot_bert_bfd")
    model = AutoModel.from_pretrained("Rostlab/prot_bert_bfd")
    #This dictionary will store proteinIDs as keys and embeddings as values
    embeddings = {}

    # Read sequences line by line. If the line starts with ">" , extract the ID; otherwise, read amino acid sequence.
    with open(input_fasta) as f:
        for line in f:
            if line.startswith(">"):
                prot_id = line.strip().lstrip(">")
            else:
                sequence = line.strip()
                #protein sequence is tokenized into input features suitable for the transformer model. return_tensors="pt" converts the input into a PyTorch tensor, the format required.
                inputs = tokenizer(sequence, return_tensors="pt")
                #This prevents computation of gradients, processes tokens and outputs hidden states for each amino acid. 
                with torch.no_grad():
                    #last_hidden_state extracts the output of the last layer in the transformer, which contains the learned representation fo reach token. mean averages the token embeddings along the sequence length dimension. squeeze().numpy() converts the pytorch tensor into a numpy array. Averaging across the entire sequence condenses the data into a single vector representing each whole protein.
                    embedding = model(**inputs).last_hidden_state.mean(1).squeeze().numpy()
                #Protein ID is key, and computed embedding (numerical vector) is the value.     
                embeddings[prot_id] = embedding

    # Save embeddings to file
    with open(output_pkl, "wb") as f:
        pickle.dump(embeddings, f)

generate_embeddings(args.input, args.output)
```

**Run the embedding script**

Embed the fasta files of the transmembrane domain containing proteins and transmembrane-less proteins. This is the most time consuming step of the pipeline, and can take >24 hours for 100k proteins. 
```
python generate_embeddings.py --input TMProteinDataset.fasta --output TMembeddings.pkl   
```

Embed the fasta files of your unknown set of fasta sequences you plan to screen for transmembrane domains. 
```
python generate_embeddings.py --input SignalPDataset.fasta --output SignalPDataset.pkl
```

### Assess embedding output
Did you get a pickle file?
```
ls -lrth *.pkl
-rw-r----- 1 rick.masonbrink proj-scinet_workshop1 4.0M Nov 20 10:25 SignalPDataset.pkl
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 316M Nov  9 08:05 TMembeddings.pkl
```

# Classification of Training Proteins

### Create labels for classification 

```
echo -e "Protein_ID\tLabel" >TMLabels.tsv
grep ">" ActualTransmembranePositive.fasta |awk '{print $0"\ttransmembrane"}' |sed 's/>//g' >>TMLabels.tsv
grep ">" ActualTransmembraneNegative.fasta |awk '{print $0"\tnon-transmembrane"}' |sed 's/>//g' >>TMLabels.tsv

```
**Copy this script to a file named: train_classifier.py**
```
import pickle
import pandas as pd
import argparse
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler


# Argument parser for shell inputs
parser = argparse.ArgumentParser(description="Train a classifier to identify transmembrane domains")
parser.add_argument('--embeddings', type=str, required=True, help="Path to embeddings pickle file")
parser.add_argument('--labels', type=str, required=True, help="Path to labels file (TSV format)")
parser.add_argument('--output', type=str, required=True, help="Path to save trained classifier pickle file")
parser.add_argument('--scaler_output', type=str, required=True, help="Path to save the scaler pickle file")
args = parser.parse_args()

# Load embeddings and labels files
with open(args.embeddings, "rb") as f:
    embeddings = pickle.load(f)

labels_df = pd.read_csv(args.labels, sep="\t")
labels_df = labels_df[labels_df["Protein_ID"].isin(embeddings.keys())]

# Loads each Protein_ID and Lable into a dataframe, converting transmembrane domain containing proteins and proteins without, zeros and ones respectively.  
X = [embeddings[prot_id] for prot_id in labels_df["Protein_ID"]]
y = labels_df["Label"].apply(lambda x: 1 if x == "transmembrane" else 0).values

# Prints the distribution of the two classes being classified to record if your datasets were balanced or not. 
num_transmembrane_proteins = sum(y)
num_non_transmembranes = len(y) - num_transmembrane_proteins
print(f"Class Distribution: {num_transmembrane_proteins} transmembrane domain-containing proteins, non_transmembrane_proteins, {num_non_transmembranes} non-transmembranes")


# Split the data into training and testing to evaluate model performance. random_state=42 ensure reproducability by fixing randomness.
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Scale the feature matrix.   This standardizes the embeddings by removin gthe mean and scaling to unit variance. The scaler computes mean and standard deviation and applies the transformation.
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train logistic regression model
model = LogisticRegression()  # Handles class imbalance
model.fit(X_train_scaled, y_train)

# Evaluate the model with Accuracy and a classification report (Precision, Recall, F1-Score, and Support). 
y_pred = model.predict(X_test_scaled)
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy}")
print("Classification Report:")
print(classification_report(y_test, y_pred))

# Save the trained classifier
with open(args.output, "wb") as f:
    pickle.dump(model, f)
print(f"Classifier saved as {args.output}")

# Save the scaler
with open(args.scaler_output, "wb") as f:
    pickle.dump(scaler, f)
print(f"Scaler saved as {args.scaler_output}")

# Print model coefficients
print("Model coefficients (weights):")
print(model.coef_)

# Print feature matrix stats
print("Feature matrix stats:")
print(f"Mean: {np.mean(X)}, Variance: {np.var(X)}")

# Confusion matrix
print("Confusion Matrix")
print("Non-transmembrane proteins correctly identified, Non-transmembrane proteins misclassified")
print("transmembrane proteins correctly identified, Transmembrane protein misclassified")
print(confusion_matrix(y_test, y_pred))
```

**Train the classifier**
```
#Here we use the embeddings we obtained from generate_embeddings.py and the labels we created to classify our sequences.
python train_classifier.py --embeddings TMembeddings.pkl --labels TMLabels.tsv --output logistic_regression_model.pkl --scaler_output scaler.pkl
```

### Evaluate your model
```
Accuracy: 0.5079077429983525
Classification Report:
              precision    recall  f1-score   support
           0       0.57      0.01      0.03      2997
           1       0.51      0.99      0.67      3073
    accuracy                           0.51      6070
   macro avg       0.54      0.50      0.35      6070
weighted avg       0.54      0.51      0.35      6070
Classifier saved as transmembraneOut.pkl
Scaler saved as scaler.pkl
Model coefficients (weights):
# Split the data into training and testing to evaluate model performance. random_state=42 ensure reproducability by fixing randomness.
[[ 0.00436481  0.00399653 -0.01072507 ...  0.00805357 -0.01355737
   0.00183267]]
Feature matrix stats:
Mean: -0.0016713324002921581, Variance: 0.009750734083354473
Confusion Matrix
Non-transmembrane proteins correctly identified, Non-transmembrane proteins misclassified
Transmembrane proteins correctly identified, Transmembrane protein misclassified
[[  41 2956]
 [  31 3042]]
```
In your stdout you should have this output from the small model. The classification report provides metrics on: <br>
*  Precision: Fraction of positive predictions that are correct.
*  Recall: Fraction of true positives that are correctly predicted.
*  F1-Score: Harmonic mean of precision and recall.
*  Support: Number of actual occurrences of each class. <br>
The zero and 1 in our table represent our classes (non-transmembrane 0) and (transmembrane 1). Ideally we would have high scores for all classes. Coefficients that are further away from zero are better, and you'd like to have a higher number for Mean and Variance. At the bottom is the 20% training data that we split being tested. 31 tramsmembrane domain containing proteins accurately identified and 3042 transmembrane containing proteins incorrectly classified, 41 nontransmembrane domain containing proteins were correct and 2956 were not correct.  

### Not good results?  Lets try a random forest model 
```
import pickle
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, precision_recall_curve, auc
from sklearn.decomposition import PCA
from imblearn.combine import SMOTETomek

# Load embeddings and labels
embeddings = pickle.load(open("TMembeddings.pkl", "rb"))
labels_df = pd.read_csv("TMLabels.tsv", sep="\t")
labels_df = labels_df[labels_df["Protein_ID"].isin(embeddings.keys())]

# Prepare input features and labels
X = [embeddings[prot_id] for prot_id in labels_df["Protein_ID"]]
y = labels_df["Label"].apply(lambda x: 1 if x == "transmembrane" else 0).values

# Convert to numpy arrays
X = np.array(X)

# Check initial class distribution
print(f"Class Distribution: {sum(y)} transmembrane, {len(y) - sum(y)} non-transmembrane")

# Apply SMOTE-Tomek
smote_tomek = SMOTETomek(random_state=42)
X_resampled, y_resampled = smote_tomek.fit_resample(X, y)
print("Class Distribution After SMOTE-Tomek:")
print(pd.Series(y_resampled).value_counts())

# Apply Edited Nearest Neighbours (ENN) if needed (otherwise, skip if no improvement)
from imblearn.under_sampling import EditedNearestNeighbours
enn = EditedNearestNeighbours()
X_resampled, y_resampled = enn.fit_resample(X_resampled, y_resampled)
print("Class Distribution After ENN:")
print(pd.Series(y_resampled).value_counts())

# Stratified train-test split to preserve class distribution
X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=42, stratify=y_resampled)

# Apply PCA for dimensionality reduction
pca = PCA(n_components=50, random_state=42)
X_train_pca = pca.fit_transform(X_train)
X_test_pca = pca.transform(X_test)
print(f"Explained variance by PCA components: {np.sum(pca.explained_variance_ratio_):.2f}")

# Train Random Forest model with class balancing
model = RandomForestClassifier(n_estimators=1000, max_depth=50, random_state=42, class_weight='balanced')
model.fit(X_train_pca, y_train)

# Make predictions and evaluate the model
y_pred = model.predict(X_test_pca)
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy}")

# Detailed classification report
print("Classification Report:")
print(classification_report(y_test, y_pred))

# Confusion matrix
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))

# Precision-Recall curve and AUC
y_proba = model.predict_proba(X_test_pca)[:, 1]  # Probability of class 1 (transmembrane)
precision, recall, _ = precision_recall_curve(y_test, y_proba)
pr_auc = auc(recall, precision)
print(f"Precision-Recall AUC: {pr_auc:.2f}")

# Save the model and PCA transformer
with open("random_forest_model.pkl", "wb") as f:
    pickle.dump(model, f)

with open("pca_transformer.pkl", "wb") as f:
    pickle.dump(pca, f)

print("Model and PCA transformer saved.")
```

### Predict transmembrane-containing proteins from your dataset of unknowns

**Copy this script to a file named: predict_transmembranes.py**
```
import pickle
import argparse
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

# Argument parser to specify the embeddings and classifier files
parser = argparse.ArgumentParser(description="Predict transmembrane domain containing proteins using a trained classifier")
parser.add_argument('--embeddings', type=str, required=True, help="Embeddings pickle file for prediction")
parser.add_argument('--model', type=str, required=True, help="Trained classifier pickle file")
parser.add_argument('--output', type=str, required=True, help="Output file to save predictions")
args = parser.parse_args()

# Load the new embeddings and the trained classifier
with open(args.embeddings, "rb") as f:
    embeddings = pickle.load(f)
with open(args.model, "rb") as f:
    model = pickle.load(f)

# Check if the scaler file exists, if not, print an error message and exit
scaler_path = "scaler.pkl"
if not os.path.exists(scaler_path):
    print(f"Error: The scaler file '{scaler_path}' was not found. Please ensure it is saved during model training.")
    exit(1)  # Exit the script if the scaler file is missing

# Load the scaler
with open(scaler_path, "rb") as f:
    scaler = pickle.load(f)

# Prepare the embeddings for prediction
# Convert the embeddings dictionary into a list of feature vectors
X_new = np.array(list(embeddings.values()))

# Scale the new embeddings (X_new) for prediction
X_new_scaled = scaler.transform(X_new)

# Predict probabilities and class labels for each protein
predictions = {}
y_pred_probs = model.predict_proba(X_new_scaled)
for prot_id, prob in zip(embeddings.keys(), y_pred_probs[:, 1]):
    # Classify based on the probability threshold
    pred = "transmembrane" if prob > 0.43 else "non-transmembrane"
    predictions[prot_id] = (pred, prob)

# Print the prediction and the probability for the transmembrane class
print("Protein\tPrediction\tTransmembrane Probability")
for prot_id, (pred, prob) in predictions.items():
    print(f"{prot_id}\t{pred}\t{prob:.4f}")

# Write predictions to the specified output file
with open(args.output, "w") as f:
    for prot_id, (label, prob) in predictions.items():
        f.write(f"{prot_id}\t{label}\t{prob:.4f}\n")
print(f"Predictions saved to {args.output}")
```

**Run the prediction**
```
python ../predict_transmembranes.py --embeddings SignalPDataset.pkl --model logistic_regression_model.pkl --output predictions.tsv
```
# Evaluate your Results
Based upon the threshold you set in your predict_transmembranes.py script, (0.43) in this case, you can determine if the protein was called a transmembrane-containing protein or non-transmembrane . This setting is directly proportional to the probability output in the third column of predictions.tsv.



SignalPDataset.fasta has 1000 proteins, 200 of which have a transmembrane domain. 
 
[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)