---
title: DeepGoPlus: Using AI to classify proteins
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


# A ProtTrans pipeline to differentiate kinases from other proteins
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
    *  Example: Classifying proteins as kinases or non-kinases using a supervised model.
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

Here I just went to Uniprot.org to find proteins that were representative kinases from UniRef (manually reviewed). Since we have been working with Arabidopsis, I have kept with that lineage by extracting only proteins from the Brassicaceae. To avoid biasing the model, we want to keep both protein types at an equal representation. 
<br>

*  Since kinases are a minute fraction of all proteins, I only extracted proteins that were in groups with at least 8 sequences with a cluster. 
```

#number of kinase sequences
grep -c ">" KinasesTax3700Uniprot.fasta
58918
```
*  For non-kinases I grabbed the cluster representative protein for each Brassicaceae cluster with at least 8 members in a cluster and Cluster name could not contain the word kinase.
```
#number of nonkinases
grep -c ">" CleanNamesNonKinasesTax3700Uniprot.fasta
20627 non-kinases
```
* Gets rid of the extra header information in the fasta files, which can affect the time to embed the proteins
```
awk '{print $1}'  uniref_NOT_name_kinase_AND_count_8_T_2024_11_11.fasta >CleanNamesNonKinasesTax3700Uniprot.fasta
awk '{print $1}'  uniref_taxonomy_id_3700_AND_name_kin_2024_11_08.fasta >CleanNamesKinasesTax3700Uniprot.fasta

cat CleanNames*fasta >ProteinDataset.fasta
```
### Create a dataset of proteins to test your trained model

```
#A small sample of kinases to test, note these kinases than what was used in training
head -n 1002 CleanNamesKinases2Tax3700Uniprot.fasta >SmallKinases.fasta
#A small sample of non-kinases, note that these non-kinases were not used in training 
head -n 19998 CleanNamesNonKinases2Tax3700Uniprot.fasta >SmallNonKinases.fasta

grep -c ">" Small*
SmallKinases2.fasta:84
SmallNonKinases2.fasta:2769

cat SmallKinases2.fasta SmallNonKinases2.fasta >SmallTrainingDataset.fasta
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

Embed the fasta files of the kinases and nonkinases. This is the most time consuming step of the pipeline, and can take >24 hours for 100k proteins. 
```
python generate_embeddings.py --input ProteinDataset.fasta --output embeddings.pkl
```

Embed the fasta files of your unknown set of fasta sequences you plan to screen for kinases. 
```
python generate_embeddings.py --input SmallTrainingDataset.fasta --output SmallTrainingembeddings.pkl
```

### Assess embedding output
Did you get a pickle file? 
```
ls -lrth *.pkl
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 12M Nov 12 12:57 SmallTrainingembeddings.pkl
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 316M Nov  9 08:05 embeddings.pkl
```

# Classification of Training Proteins

### Create labels for classification 
```
#creates header
echo -e "Protein_ID\tLabel" >labels.tsv

#grabs the fasta header and appends a print of fasta header name tab kinase to labels.tsv
grep ">" KinasesTax3700Uniprot.fasta |awk '{print $1"\tkinase"}' |sed 's/>//g' >>labels.tsv

#grabs the fasta header and appends 'fasta_header_name tab non-kinase' to labels.tsv
grep ">" NonKinasesTax3700Uniprot.fasta |awk '{print $1"\tnon-kinase"}' |sed 's/>//g' >>labels.tsv
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
parser = argparse.ArgumentParser(description="Train a classifier to identify kinases")
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

# Loads each Protein_ID and Lable into a dataframe, converting kinase and non-kinase to zeros and ones. 
X = [embeddings[prot_id] for prot_id in labels_df["Protein_ID"]]
y = labels_df["Label"].apply(lambda x: 1 if x == "kinase" else 0).values

# Prints the distribution of the two classes being classified to record if your datasets were balanced or not. 
num_kinases = sum(y)
num_non_kinases = len(y) - num_kinases
print(f"Class Distribution: {num_kinases} kinases, {num_non_kinases} non-kinases")


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
print("Non-kinase correctly identified, Non-kinase misclassified")
print("Kinase correctly identified, Kinase misclassified")
print(confusion_matrix(y_test, y_pred))
```

**Train the classifier**
```
#Here we use the embeddings we obtained from generate_embeddings.py and the labels we created to classify our sequences.
python train_classifier.py --embeddings embeddings4PredictionDataset.pkl --labels labels4PredictionDataset.tsv --output logistic_regression_model.pkl --scaler_output scaler.pkl
```

### Evaluate your model
```
              precision    recall  f1-score   support

           0       0.98      1.00      0.99       558
           0       0.98      1.00      0.99       558
           1       0.00      0.00      0.00        13

    accuracy                           0.98       571
   macro avg       0.49      0.50      0.49       571
weighted avg       0.95      0.98      0.97       571

Classifier saved as logistic_regression_modelSmall.pkl
Scaler saved as scalerSmall.pkl
Model coefficients (weights):
[[-0.01420174 -0.00063592 -0.00290767 ... -0.00020602  0.00301244
  -0.00229889]]
Feature matrix stats:
Mean: -0.0016713655786588788, Variance: 0.009751115925610065
Confusion Matrix
Non-kinase correctly identified, Non-kinase misclassified
Kinase correctly identified, Kinase misclassified
[[558   0]
 [ 13   0]]

In your stdout you should have this output from the small model. We are looking for high F1 scores, coefficients that are further away from zero are better, and you'd like to have a higher number for Mean and Variance.  At the bottom is the 20% training data that we split being tested.  558 nonkinases were correct and 13 kinases were correct, without any misclassifications. 
```

### Predict kinases from your dataset of unknowns

**Copy this script to a file named: predict_kinases.py**
```
import pickle
import argparse
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

# Argument parser to specify the embeddings and classifier files
parser = argparse.ArgumentParser(description="Predict kinases using a trained classifier")
parser.add_argument('--embeddings', type=str, required=True, help="Embeddings pickle file for prediction")
parser.add_argument('--model', type=str, required=True, help="Trained classifier pickle file")
parser.add_argument('--output', type=str, required=True, help="Output file to save predictions")
args = parser.parse_args()

# Load the new embeddings 
with open(args.embeddings, "rb") as f:
    embeddings = pickle.load(f)

# load the trained classifier
with open(args.model, "rb") as f:
    model = pickle.load(f)

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

# classify your prediction based on the probability that it is a kinase
for prot_id, prob in zip(embeddings.keys(), y_pred_probs[:, 1]):
    # Classify based on the probability threshold
    pred = "kinase" if prob > 0.6 else "non-kinase"
    predictions[prot_id] = (pred, prob)

# Print the prediction and the probability for the kinase class
print("Protein\tPrediction\tKinase Probability")
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
python predict_kinases.py --embeddings SmallTrainingembeddings.pkl --model logistic_regression_model.pkl --output predictions.tsv
```
# Evaluate your Results

[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)