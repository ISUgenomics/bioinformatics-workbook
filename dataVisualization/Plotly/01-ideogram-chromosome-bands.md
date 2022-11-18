---
title: "Ideogram: display chromosome bands"
layout: single
author: Aleksandra Badaczewska
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

{% include toc %}

# Overview

An Ideogram, available in the <a href="https://mail.google.com/mail/u/0/?zx=r2gbow53bckx#inbox" target="_blank">dash-bio</a> set of visualization approaches, is a **diagram composed of multiple bars** displayed vertically or horizontally. The range of positions/indexes for a measured feature determine the relative length of the bar. A set of discrete **colors defines the type or intensity of a feature** along the graph.

![ideogram](../assets/images/ideogram.png)

# Applications

In ***Bioinformatics*** ideogram graph is used to **visualize the positions of genes** or microRNA along the chromosome. Easily distinguishable regions are called cytogenetic bands. The **pattern of bands** makes each chromosome unique. The ideogram shows the set of chromosomes of varying lengths, every split into two arms divided by a centromere.

![ideogram](../assets/images/ideogram_bioinformatics.png)

Using an ideogram, you can also ***visualize the results of electrophoresis***. It is an experimental technique where the biomolecules (peptides, proteins, nucleic acids) are segregated **based on their size and electrical charge**. Depending on these properties, the **molecules migrate** through the matrix under electric current at different rates, **leaving the characteristic pattern**. Analyzing the **location and size of bands** allows for identifying the type and number of molecules.

![ideogram](../assets/images/ideogram_electrophoresis.png)

Although the mentioned are the most common applications of the ideogram, you can also use this chart type to **visualize other data**. Typically, it is a good choice for visualizing **irregular ranges along an ordered direction**.
The **direction** can be an observation index, position, time, or any unitary variable. The difference between the highest and smallest value will **determine the length of the bar**.
Ranges or **bands** along the bar **visualize the pattern** of the analyzed feature. Different **colors** may correspond to labels, categories, or **feature types**. *For example, they can depict the species of grains crop that has given the highest yields over the decades in the last century.* But also, different **colors** can refer to the **intensity of a feature**. *For example, they can show how the amount of corn produced has changed over the decades in the last century.*

**Table.** The yields *[in 1000 MT]* of grain crops in the United States, aggregated over decades using <a href="https://github.com/ISUgenomics/data_wrangling" target="_blank">ISUgenomics/data_wrangling </a>/ <a href="https://github.com/ISUgenomics/data_wrangling/tree/main/bin_data" target="_blank">bin_data  ⤴</a> mini python application *(usage explained in this tutorial)*.

```
     year barley    corn   oats  rice  rye  wheat
-------------------------------------------------
1960-1969  87962 1038452 136986 25047 7645 359652
1970-1979  88053 1513440  98472 34732 6345 496421
1980-1989 106060 1817012  63564 46367 5089 637308
1990-1999  83818 2192185  30834 56466 2626 644425
2000-2009  52666 2795878  16296 66408 1899 570611
2010-2019  40055 3425655   9309 64921 2265 557878
```

<details><summary> The raw data taken from the source of: <a href="https://www.indexmundi.com/agriculture/?country=us&commodity=milled-rice&graph=production" target="_blank">Index Mundi: Agriculture  ⤴</a> <i>[here: example <b>rice</b> data]</i> </summary>

<br>
Year    Production    Unit_of_Measure Growth_Rate <br>
1960    1756    (1000 MT)       NA      <br>
1961    1763    (1000 MT)       0.40%   <br>
1962    2133    (1000 MT)       20.99%  <br>
1963    2295    (1000 MT)       7.59%   <br>
1964    2386    (1000 MT)       3.97%   <br>
1965    2497    (1000 MT)       4.65%   <br>
1966    2805    (1000 MT)       12.33%  <br>
1967    2950    (1000 MT)       5.17%   <br>
1968    3459    (1000 MT)       17.25%  <br>
1969    3003    (1000 MT)       -13.18% <br>
1970    2796    (1000 MT)       -6.89%  <br>
1971    2838    (1000 MT)       1.50%   <br>
1972    2828    (1000 MT)       -0.35%  <br>
1973    3034    (1000 MT)       7.28%   <br>
1974    3667    (1000 MT)       20.86%  <br>
1975    4099    (1000 MT)       11.78%  <br>
1976    3781    (1000 MT)       -7.76%  <br>
1977    3120    (1000 MT)       -17.48% <br>
1978    4271    (1000 MT)       36.89%  <br>
1979    4298    (1000 MT)       0.63%   <br>
1980    4810    (1000 MT)       11.91%  <br>
1981    5979    (1000 MT)       24.30%  <br>
1982    4960    (1000 MT)       -17.04% <br>
1983    3215    (1000 MT)       -35.18% <br>
1984    4382    (1000 MT)       36.30%  <br>
1985    4332    (1000 MT)       -1.14%  <br>
1986    4307    (1000 MT)       -0.58%  <br>
1987    4109    (1000 MT)       -4.60%  <br>
1988    5186    (1000 MT)       26.21%  <br>
1989    5087    (1000 MT)       -1.91%  <br>
1990    5098    (1000 MT)       0.22%   <br>
1991    5096    (1000 MT)       -0.04%  <br>
1992    5704    (1000 MT)       11.93%  <br>
1993    5053    (1000 MT)       -11.41% <br>
1994    6384    (1000 MT)       26.34%  <br>
1995    5628    (1000 MT)       -11.84% <br>
1996    5453    (1000 MT)       -3.11%  <br>
1997    5750    (1000 MT)       5.45%   <br>
1998    5798    (1000 MT)       0.83%   <br>
1999    6502    (1000 MT)       12.14%  <br>
2000    5941    (1000 MT)       -8.63%  <br>
2001    6714    (1000 MT)       13.01%  <br>
2002    6536    (1000 MT)       -2.65%  <br>
2003    6420    (1000 MT)       -1.77%  <br>
2004    7462    (1000 MT)       16.23%  <br>
2005    7101    (1000 MT)       -4.84%  <br>
2006    6267    (1000 MT)       -11.74% <br>
2007    6288    (1000 MT)       0.34%   <br>
2008    6546    (1000 MT)       4.10%   <br>
2009    7133    (1000 MT)       8.97%   <br>
2010    7593    (1000 MT)       6.45%   <br>
2011    5866    (1000 MT)       -22.74% <br>
2012    6348    (1000 MT)       8.22%   <br>
2013    6117    (1000 MT)       -3.64%  <br>
2014    7106    (1000 MT)       16.17%  <br>
2015    6131    (1000 MT)       -13.72% <br>
2016    7117    (1000 MT)       16.08%  <br>
2017    5659    (1000 MT)       -20.49% <br>
2018    7107    (1000 MT)       25.59%  <br>
2019    5877    (1000 MT)       -17.31% <br>
2020    7224    (1000 MT)       22.92%  <br>
2021    6090    (1000 MT)       -15.70% <br>
2022    5589    (1000 MT)       -8.23%

</details><br>


![Grain crops](../assets/images/ideogram-crops.png)

<br><br>
# Case study

Let's assume the corn yield **data at daily frequency was collected** (*^365 times a year at most*) for a hundred years with continuing indexing of days. Now, we are interested in **identifying general periods of shortages**. A change in the data structure is required before visualization to get an informative output. To begin with, we need to **aggregate the data to make a more coarse-grained unit of time**.

![Data aggregation](../assets/images/ideogram-data-aggregation.png)

***How will we aggregate the data?*** <br>
Most simply by summing or averaging the data from the selected period.

***How to optimize the length of the period?*** <br> Ideally, the way to highlight the significant level of feature variability.

1) In the general scenario, we could aggregate the data over every 365 input rows, slicing it into 100 annual periods. <br>
2) In another variant, we could ask for ten periods only and estimate the number of rows merged into the data slice.

**However, let's suppose that data from some random days over the years are missing**. <br>
In this case, dividing the data into **chunks based on a fixed number of rows or slices will lose the reference to time**. (*In the first scenario, the last period will be much smaller, while in the second scenario, all periods will be equal but still smaller than 365 days, and some days can drop into the wrong year*). <br>

![Data slicing](../assets/images/ideogram-data-slices.png)

So, to solve this issue, we should **create data slices based on the increment** of the day counter, where days from the first year to the 100th year were indexed continually, including days of missing data. <br>
Thus, we will take all rows whose indexes match the increment range for every slice. The value increment should be 365 to create annual periods. *(So, 1-365 are days indexes for the first year, 366 - 730 for the second year, etc.).* <br>

![](../assets/images/ideogram-data-slices-good.png)

**Note that in the slices where missing yields occurred, the count of the rows (days) is smaller but known.** Thus, since the number of observations in periods varies, **data aggregation by averaging** over the data slice seems a more robust solution.

<br><br>
# Hands-on tutorial

In this practical tutorial we will use the bioinformatic data representing genetic features of species X detected in the chromosomes of newly discovered organism Y.

## **Raw data**

The data collected is a matrix indexed by a chromosome ID (*first column*) and positions along it (*second column*). The following numerical columns correspond to several (8) individuals of Y species. They contain an alignment depth given per position in the chromosome.

*Preview of `raw_data.txt`* file
<a href="https://downgit.github.io/#/home?url=https://github.com/ISUgenomics/data_wrangling/blob/main/bin_data/bin_data.py" target="_blank"><input type="button" value="Download .zip ⤵" style="background-color: #e7e7e7; color: black;" /></a>
<a href="https://raw.githubusercontent.com/ISUgenomics/bioinformatics-workbook/master/dataVisualization/Plotly/01-ideogram/raw_data.txt" target="_blank"><input type="button" value="Open raw ⤴" style="background-color: #e7e7e7; color: black;" /></a>
<a id="raw-url" href="https://github.com/ISUgenomics/bioinformatics-workbook/tree/master/dataVisualization/Plotly/01-ideogram" target="_blank"><input type="button" value="Open at GitHub ⤴" style="background-color: #e7e7e7; color: black;" /></a>

```
HiC_scaffold_1  982     0       0       0       0       0       1       0       0
HiC_scaffold_1  983     0       0       0       0       0       1       0       0
HiC_scaffold_1  984     0       0       0       0       0       1       0       0
HiC_scaffold_1  985     0       0       0       0       0       1       0       0
HiC_scaffold_1  986     0       0       0       0       0       1       0       0
...
```

Ten chromosomes are of interest to us in total, each of the lengths (*number of positions*) specified in the second column.

*Preview of `scaffold_list`* file
<a href="https://raw.githubusercontent.com/ISUgenomics/bioinformatics-workbook/master/dataVisualization/Plotly/01-ideogram/scaffold_list" target="_blank"><input type="button" value="Open raw ⤴" style="background-color: #e7e7e7; color: black;" /></a>
<a id="raw-url" href="https://github.com/ISUgenomics/bioinformatics-workbook/blob/master/dataVisualization/Plotly/01-ideogram/scaffold_list" target="_blank"><input type="button" value="Open at GitHub ⤴" style="background-color: #e7e7e7; color: black;" /></a>

```
HiC_scaffold_1 73700
HiC_scaffold_2 46643
HiC_scaffold_3 39042
HiC_scaffold_4 17063
HiC_scaffold_5 47495
HiC_scaffold_6 44324
HiC_scaffold_7 114593
HiC_scaffold_9 13668
HiC_scaffold_8 22968
HiC_scaffold_10 41740
```

Assuming the raw file comes contaminated with results for other scaffolds, we will keep only the rows that match the chromosome list for the Y organism.
```
for i in `cat scaffold_list | awk '{print $1}'`; do grep $i"\t" < raw_data.txt >> input_data.txt; done
```

The resulting file `input_data.txt` will be used as a direct input for the data aggregation step.

## **Bin Data app**

To process the `input_data.txt` file we use a ready-made Python application <b><a href="https://github.com/ISUgenomics/data_wrangling/bin_data" target="_blank">bin_data ⤴</a></b>, available in the <a href="https://github.com/ISUgenomics/" target="_blank"> ISUgenomics </a>/<a href="https://github.com/ISUgenomics/data_wrangling" target="_blank">data_wrangling ⤴</a> repository.

Get `bin_data.py` script directly
<a href="https://downgit.github.io/#/home?url=https://github.com/ISUgenomics/data_wrangling/blob/main/bin_data/app/bin_data.py" target="_blank"><input type="button" value="Download .zip ⤵" style="background-color: #e7e7e7; color: black;" /></a>
<a href="https://raw.githubusercontent.com/ISUgenomics/data_wrangling/main/bin_data/app/bin_data.py" target="_blank"><input type="button" value="Open raw ⤴" style="background-color: #e7e7e7; color: black;" /></a>
<a id="raw-url" href="https://github.com/ISUgenomics/data_wrangling/blob/main/bin_data/app/bin_data.py" target="_blank"><input type="button" value="Open at GitHub ⤴" style="background-color: #e7e7e7; color: black;" /></a>

You can learn more about the <a href="https://github.com/ISUgenomics/data_wrangling/tree/main/bin_data" target="_blank">bin_data.py ⤴</a> app from the comprehensive tutorial <a href="https://datascience.101workbook.org/07-DataParsing/03-DATA-WRANGLING-APPS/02-slice-or-bin-data-py" target="_blank">Aggregate data over slicing variations ⤴</a> published in the <a href="https://datascience.101workbook" target="_blank">Data Science Workbook ⤴</a>. You will find there:
* the application <a href="https://datascience.101workbook.org/07-DataParsing/03-DATA-WRANGLING-APPS/02-slice-or-bin-data-py#app-overview" target="_blank">Overview ⤴</a>
* description of the <a href="https://datascience.101workbook.org/07-DataParsing/03-DATA-WRANGLING-APPS/02-slice-or-bin-data-py#algorithm/" target="_blank">Algorithm ⤴</a>
* documentation, including <a href="https://datascience.101workbook.org/07-DataParsing/03-DATA-WRANGLING-APPS/02-slice-or-bin-data-py#requirements" target="_blank">Requirements ⤴</a> and <a href="https://datascience.101workbook.org/07-DataParsing/03-DATA-WRANGLING-APPS/02-slice-or-bin-data-py#options" target="_blank">Options ⤴</a>
* examples of generic <a href="https://datascience.101workbook.org/07-DataParsing/03-DATA-WRANGLING-APPS/02-slice-or-bin-data-py#usage-generic" target="_blank">Usage ⤴</a>
* <a href="https://datascience.101workbook.org/07-DataParsing/03-DATA-WRANGLING-APPS/02-slice-or-bin-data-py#hands-on-tutorial" target="_blank">Hands-on tutorial ⤴</a>, including 6 detailed case studies

### ❖ *General Approach*

Before we aggregate data over data slices, first, we need to split it into **chromosome-based data chunks**. To do so, we can use the column with chromosome IDs, referred further to as `label-col`. Then, each data chunk will split into **custom-size slices using the increment** of values of the column storing positions, referred further to as `ranges-col`. Finally, numerical `data-cols` (columns of different traits) will aggregate for every slice by calculating the `sum` or `mean`.

![Algorithm](../assets/images/bin_data.png) <br>
<i>The figure shows the main steps of the <code>bin_data</code> algorithm.</i><br>

The optimal data structure requires:<br>
<b>L</b> - <code>label-col</code>, a column of labels, <b>HERE:</b> chromosome IDs<br>
<b>R</b> - <code>ranges-col</code>, numerical [int, float] column of ordered data characteristic, <b>HERE:</b> positions along the chromosome<br>
<b>numD</b> - <code>data-cols</code>, any number of numerical columns that will be aggregated, <b>HERE:</b> traits 1-8.<br>

![Data structure](../assets/images/bin_data_structure.png) <br>
<i>The figure shows the data structure of `input_data.txt`. The <b>L</b> column stores chromosomes IDs and is used to create label-based data chunks. The <b>R</b> column stores positions along the chromosome and is used to create slicing ranges. The <b>numD</b> columns are numerical data that will be aggreagted over slices.</i>

<span style="color: #ff3870;font-weight: 500;">This step is implemented in the <b>bin_data.py</b> Python app.</span>

### ❖ *Aggregate observations over value increment*

When you look closely at the data, you notice that **some positions** along the chromosomes **are missing** in the matrix.

![Missing data](../assets/images/bin_data_missing.png)

Because of that, aggregating data over <u>a constant number of rows</u> `[-t 'step']` or splitting data into <u>a fixed number of slices</u> `[-t 'bin']` will lose the reference to exact positions along the chromosome. Thus for such a dataset, we **aggregate observations over value increment** `[-t 'value']`, where the value is the positions column. <br>
In this case, we request to aggregate data with increment X `[-n X]` of position values. *For example, with increment X=1000, the first data slice is created from position 1 to 1000, the second from position 1001 to 2000, etc.* To detect the most **optimal slicing increment**, we should consider a few variants of X. Note the average length of the chromosomes in the input file (*significantly truncated*) is several tens of thousand, so a good choice is setting X to 100, 1k, and 10k (*would be millions for a full-length chromosome*). <br>
Considering data for some indexes are missing, the number of input rows **(counts) may vary among slices**. Since the lack of data is random, we can not be sure if these positions would be significant in detecting the trait. Therefore, it makes more sense to **average data over slices** `[-c 'ave']` instead of summing it `[-c 'sum']`. To make sure numerical results round to a meaningful number, let us keep three decimal places `[-d 3]`.

***Run `bin_data.py` app in the command line:***

Get `bin_data.py` script:
<a href="https://downgit.github.io/#/home?url=https://github.com/ISUgenomics/data_wrangling/blob/main/bin_data/app/bin_data.py" target="_blank"><input type="button" value="Download .zip ⤵" style="background-color: #e7e7e7; color: black;" /></a>
<a href="https://raw.githubusercontent.com/ISUgenomics/data_wrangling/main/bin_data/app/bin_data.py" target="_blank"><input type="button" value="Open raw ⤴" style="background-color: #e7e7e7; color: black;" /></a>
<a id="raw-url" href="https://github.com/ISUgenomics/data_wrangling/blob/main/bin_data/app/bin_data.py" target="_blank"><input type="button" value="Open at GitHub ⤴" style="background-color: #e7e7e7; color: black;" /></a>
on your local machine.

Run the program following the command:
```
python3 bin_data.py -i input_data.txt -l 0 -r 1 -t 'value' -n 100 -c 'ave' -d 3 -o 'output-value_ave_100'
```

In addition to the aggregated data output (*output-value_ave_100.csv*), the application will generate a `CHUNKS` directory where it writes **chromosome-based chunks of raw data**. You can use these files as input in the repetitions of the analysis, such as in the X-increment optimization. That will definitely speed up the processing of Big Data inputs.

*Preview or get the `CHUNKS`* folder
<a href="https://downgit.github.io/#/home?url=https://github.com/ISUgenomics/bioinformatics-workbook/tree/master/dataVisualization/Plotly/01-ideogram/CHUNKS" target="_blank"><input type="button" value="Download .zip ⤵" style="background-color: #e7e7e7; color: black;" /></a>
<a id="raw-url" href="https://github.com/ISUgenomics/bioinformatics-workbook/tree/master/dataVisualization/Plotly/01-ideogram/CHUNKS" target="_blank"><input type="button" value="Open at GitHub ⤴" style="background-color: #e7e7e7; color: black;" /></a>

```
for X in 1000 10000
do
    python3 bin_data.py -i ./CHUNKS -l 0 -r 1 -t 'value' -n $X -c 'ave' -d 3 -o 'output-value_ave_$X'
done
```

*Preview of the `output-value_ave_10000` with data aggregated over slices cut with X=10k units in the 'position' column. The trait's values were averaged for rows in every data slice, and the number of rows is given in the 'counts' column.*

```
         label,   position,val-2,val-3,val-4,val-5,val-6,val-7,val-8,val-9,count
HiC_scaffold_1,  982-10982,0.591,0.567,0.048,0.000,0.380,0.822,0.074,0.359,1373
HiC_scaffold_1,10982-20982,0.803,0.642,0.238,0.326,0.895,1.057,0.000,0.456,2025
HiC_scaffold_1,20982-30982,0.107,0.000,0.000,0.052,0.507,0.378,0.000,0.085,672
HiC_scaffold_1,30982-40982,0.471,0.000,0.417,0.000,0.206,0.007,0.000,0.000,427
HiC_scaffold_1,40982-50982,0.865,0.697,0.152,0.228,0.410,1.276,0.151,1.030,1279
HiC_scaffold_1,50982-60982,0.313,0.217,0.332,0.082,0.362,0.497,0.034,0.033,2042
HiC_scaffold_1,60982-70982,1.580,1.008,0.204,0.288,0.776,1.593,0.455,0.195,1768
```


## **Visualize using ideogram**

Once the data is aggregated to different coarseness, it is time to **adjust the format for visualization on the ideogram.** <br>
An ideogram is a diagram in which we have only one dimension for a numerical variable. This is the index or position along the bar. <br>
***So how can you display aggregated data for individual traits?*** <br>
We need to convert the numerical values into colors used for regions bounded by the position ranges (referred further as **bands**). <br>
***How to do it?*** <br>
First you need to choose a scale composed of discrete colors. Then, each color needs to be assigned a range of values that correspond to it. <br>
***How can you map values to colors?*** <br>
First, you need to find the minimum, maximum, mean, and standard deviation values for the trait's aggregated data. Assign values around the mean for the color corresponding to the center of the scale. Then, gradually adjust values between minimum and mean for the lower part of the color scale and values between mean and maximum for the upper part of the color scale. The number of value ranges should correspond to the number of colors used. You can also depend on the standard deviation for the length of the value range for a given color.

![Values to colors mapping](../assets/images/ideogram_convert_data.png) <br>
*The figure shows the values to colors mapping.*

**Luckily, you do NOT have to implement the algorithm that will do it yourself.** <br>
<span style="color: #ff3870;font-weight: 500;">Use a ready-made <a href="https://github.com/ISUgenomics/data_wrangling/blob/main/bin_data/convert_for_ideogram.py" target="_blank">convert_for_ideogram.py ⤴</a> application in Python</span>, moving to subsection: <br>
<b>A.</b> [dash-bio ideogram](https://bioinformaticsworkbook.org/dataVisualization/Plotly/01-ideogram-chromosome-bands.html#a-use-dash-bio-variant-js) or <br>
<b>B.</b> [customizable plotly ideogram](https://bioinformaticsworkbook.org/dataVisualization/Plotly/01-ideogram-chromosome-bands.html#b-use-plotly-variant) <br>
to follow detailed instructions for converting data for visualization on an ideogram.

### ***A. Use dash-bio variant***

An <a href="https://dash.plotly.com/dash-bio/ideogram" target="_blank">Ideogram ⤴</a> visualization is available among components in the <a href="https://dash.plotly.com/dash-bio" target="_blank">Dash Bio ⤴</a> module of the <a href="https://dash.plotly.com/introduction" target="_blank">Dash ⤴</a> framework. The <a href="https://github.com/eweitz/ideogram" target="_blank">Ideogram ⤴</a> was originally developed in JavaScript by <a href="https://github.com/eweitz" target="_blank">Eric Weitz ⤴</a> to provide chromosome visualization for the web. Further, it was encapsulated into the Dash framework for easy use by importing from the dash-bio library.

An ideogram imported from the dash-bio library makes it easy to visualize the pattern of bands on chromosomes. The application allows **visualizing chromosomes for standard organisms** available in an online <a href="https://unpkg.com/browse/ideogram@1.37.0/dist/data/" target="_blank">UNPKG database ⤴</a>. For these organisms, you can use your own annotation file. However, the approach is browser-based, and both bands and annotations **data are loaded exclusively through the URL**. That is a <u>huge limitation</u> for individual users who would like <u>to load input directly from the file system</u> on the local machine. To visualize your own data for a customized organism, you need to create an online database (*using GitHub Pages as the simplest solution*) or add files to the existing one to load data from it via URL. Another complication is the **strictly defined format** of the bands and annotations data files. Also, there are very **limited options for coloring bands**, which provides only 7 built-in gray shades (including white and black) and 2 additional colors (pink, and purple). Still, this variant of the ideogram allows for a nice interactive transition between **horizontal or vertical orientation** and **adjusting the dimensions** of the bars.

![Ideogram JS](../assets/images/ideogram_overview.gif)

To facilitate visualizations using the JavaScript variant of the Ideogram, I built a web application using the Dash-Bio ideogram component and Dash framework. The app runs within Jupyter Notebook, and a single click opens the Ideogram interface in the separate browser tab. Then, you can select inputs and adjust visualization options using Dash widgets without the necessity of any changes in the source code. <br>
Following the steps in this section, you will learn how to: <br>
(1) convert inputs into the required format <br>
(2) upload data into the online database <br>
(3) run the app in the JupyterLab <br>
(4) adjust visualization using the Ideogram

<span style="color: #ff3870;font-weight: 500;">If you want to load customized files directly from the local file system, go to section B. <a href="https://bioinformaticsworkbook.org/dataVisualization/Plotly/01-ideogram-chromosome-bands.html#b-use-plotly-variant" >Customizable Plotly-based Ideogram</a>.</span>

### ❖ *Convert data structure*

**Chromosome bands data**

The required input format for chromosome bands data is a JSON string composed of a dictionary with a `"chrBands"` keyword for which the value is a list of bands for all chromosomes.

```
{"chrBands" : []}
```
Every band is defined as a string of 8 fields separated by a single space. Subsequent strings for bands are separated by a comma.

```
 1 2    3   4        5   6        7      8
------------------------------------------
"1 p p1-0 982 10000982 982 10000982 gpos33",
"1 q q1-0 235 10000235 235 10000235 gpos50"
```

**1)** The first column is a chromosome ID, typically a number or letter (*e.g., for a human 1-22, X, Y*). <br>
**2)** The second column is a chromosome arm, `p` or `q` for shorter and longer arm, respectively. <br>
**3)** The third field is a custom name for a band, typically a gene name. <br>
**4-7)** The fourth and fifth columns contains `from` and `to` range of chromosome positions. Usually the same range is repeated in the sixth and seventh columns. <br>
**8)** The last field contains alias for the color of the band.

**Convert bands data for ideogram**

Use <a href="https://github.com/ISUgenomics/data_wrangling/blob/main/bin_data/convert_for_ideogram.py" target="_blank">convert_for_ideogram.py ⤴</a> Python mini app to assign colors for value ranges. Implemented options make the application flexible to user needs. You can learn more from the documentation available in the <a href="https://github.com/ISUgenomics/data_wrangling" target="_blank">data_wrangling</a>/<a href="https://" target="_blank">assign_colors</a> repo and the <a href="https://" target="_blank">X</a> tutorial in the <a href="https://" target="_blank">Data Science</a> workbook.

*available options of `convert_for_ideogram.py`:*

```
  -h,       --help                       show this help message and exit
  -i input, --data-source input [string] input multi-col file
  -l label, --labels-col label     [int] index of column with labels
  -r range, --ranges-col range     [int] index of column with ranges
  -a arm,   --arms-col arm         [int] index of column with chromosome arms annotation
  -b band,  --bands-col band       [int] index of column with bands annotation (names)
  -v vals,  --values-col vals      [int] list of indices of numerical columns to have color assigned
```

The `-i`, `-l`, and `-r` options are required, so you have to specify the input file, the index of a column with chromosome labels, and the index of a column with the band's position ranges. *Note numbering in Python starts from 0.* If you do NOT specify the `arm` column, the 'p' value will be assigned for all bands, and so all chromosomes will consist of a single arm. If you do NOT specify the `band` column, names of bands will be created automatically following the syntax `arm+chromosomeID-{next_int}`. If you do NOT provide the 'vals' list with indexes of traits columns, then by default, all numerical columns will be mapped to colors.

In our dataset, we do NOT have the information about arms and bands names. We also want to make value-to-color mapping for all traits. So, let it all be default for the purpose of data format requirements. We only specify the `label-col` as 0 and `ranges-col` as 1.

```
python3 convert_for_ideogram.py -i output-value_ave_10k.csv -l 0 -r 1
```

By default, the algorithm uses the `mean` value of a given trait to specify the ranges matching the built-in colors. So, the gpos50 corresponds to the mean, and gpos100 corresponds to the 2x mean. Pink and purple highlight bands with the highest intensity of a trait. <i>To learn more about other variants of value-to-color mapping, see the <a href="https://" target="_blank">assign_color.py ⤴</a> documentation.</i>

![Ideogram JS colrscale](../assets/images/ideogram_js_colorscale.png)

As the output, you will get separate JSON files for all traits (*8 in this example*).

*Preview of the resulting `data-val2-10M.json` file:*

```
{"chrBands":
  [
    "1 p p1-0 982 10000982 982 10000982 gpos33",
    "1 p p1-1 10000982 20000982 10000982 20000982 gpos50",
    "1 p p1-2 20000982 30000982 20000982 30000982 gpos33",
    ...
    "2 p p2-9 285 10000285 285 10000285 gpos66",
    "2 p p2-10 10000285 20000285 10000285 20000285 gpos33",
    "2 p p2-11 20000285 30000285 20000285 30000285 gpos25"
    ...
  ]
}
```

### ❖ *Upload data into the database*

The next step is to upload JSON files into the online database, to access them via URL. You can add files to any website/server that provides online visibility. For individual users, the simplest solution is to create <a href="https://docs.github.com/en/pages/getting-started-with-github-pages/about-github-pages" target="_blank">GitHub Pages  ⤴</a> from your personal account and create indexed file struture, that will display content online.

Once you have personal <a href="https://docs.github.com/en/get-started/signing-up-for-github/signing-up-for-a-new-github-account" target="_blank">account on GitHub ⤴</a>, please follow the instructions in the <a href="https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site" target="_blank">Creating a GitHub Pages site ⤴</a> tutorial, to host a website with your database.

<div style="background: mistyrose; padding: 15px; margin-bottom: 20px;">
<span style="font-weight:800;">WARNING:</span>
<br><span style="font-style:italic;">
<i><b>Note the website will be public, so the content (including your files) will be visible to anyone.</b> If you want to visualize data directly from your file system using local server (incognito), go to section B. <a href="https://" >Customizable Plotly-based Ideogram</a>.</i>
</span>
</div><br>

Once your website is hosted online, `git clone` the X repository on your local machine, copy the ideogram folder into yor **user.github.io** repository, and `git push` it back on the GitHub server. Follow instructions below, to do that step-by-step.

0) If you haven't already done so, clone **your** `user.github.io` repository.

![Clone repo](../assets/images/clone_repo.png) <br>
*The Figure shows how to copy the URL link of a selected GitHub repository.*

In the terminal window at selected location in the file system on your local machine, type:
```
git clone {copied URL}
```

1) Clone ideogram_db repository on your local machine.

```
git clone https://github.com/ISUgenomics/ideogram_db.git
```

The `ls` command should display both new repositories on your current path.

2) Enter the `ideogram_db` repo, remove `.git` folder, and move the remaining content into your personal repository

```
cd ideogram_db
rm -r .git
cp -r * /path/to_your_personal_repo.github.io
cd ../
rm ideogram_db
```

3) Navigate to your `.github.io` folder and make `git push`

```
cd path/{your}.github.io
git status
git add *
git commit -m 'added basic data structure of ideogram database'
git push
```

Visit your GitHub Pages online at `https://{account}.github.io` and refresh the website. Compare the appearance of the page with the expected result.

![Ideogram database](../assets/images/ideogram_database.png)

**Copy custom data into the database**

In the terminal window, copy your customized JSON data on the `ideogram/database/bands` path in your `{user}.github.io`:

```
cp path/to/JSON/data/*.json {user}.github.io/ideogram/database/bands/
```

Then, navigate to `{user}.github.io/ideogram/` and execute `run_me_before_commit.sh` bash script. This step will automatically update the index of files in the database and rename them to follow the syntax `custom-{int}` required by Ideogram JS module. *Don't worry, you will see the original names of your files in the application interface. The list of matching pairs is in the `list` file on the `.../ideogram/database/bands` path in your `{user}.github.io`.*

```
cd {user}.github.io/ideogram/
. ./run_me_before_commit.sh
```

The final step is to commit the updated database:

```
git status
git add {copy list of added and modified files}
git push
```

You can review the list of files hosted in your database online at `https://{account}.github.io/ideogram/database/bands/`.

### ❖ *Open ideogram in JupyterLab*

**Requirements**

* `jupyterlab` - web-based programming environment
* `dash`       - library of widgets for web-based applications
* `dash_bio`   - library of interactive graphs for biology tasks

The **Dash-Bio graphing library** has incorporated the original Ideogram module, written in JavaScript. That made it possible to call the diagram within web applications that uses **efficient Python data wrangling** in the back-end and **user-friendly Dash widgets** (options components) in the interactive front-end interface. My application is developed within the **JupyterLab environment**, making it robust, transferable, and easy to use daily across operating systems and web browsers.

<div style="background: #cff4fc; padding: 15px;">
<span style="font-weight:800;">PRO TIP:</span>
<br><span style="font-style:italic;"> If you are not familiar with JupyterLab yet, please follow the tutorial <a href="https://datascience.101workbook.org/04-DevelopmentEnvironment/02A-jupyter-basics" target="_blank">Jupyter: Web-Based Programming Interface  ⤴</a> available in the <a href="https://datascience.101workbook.org" target="_blank">Data Science Workbook  ⤴</a></span>.<br><br>
You can <b>install jupyterlab</b> using <b>pip</b>: <br>
$ <code>pip install jypyterlab</code> <br><br>
Then, you can launch the web-based interface from the command line: <br>
$ <code>jupyter lab</code> <br><br>
That will open your Jupyter session in a web browser on localhost with a default URL: <i><b>http://localhost:8889/lab</b></i>.
</div><br>

You can also install the remaining requirements (*dash, dash_bio*) with `pip` directly in the terminal, making them available to the entire system. However, using the **Conda environment manager** will be a neater solution.

<div style="background: #cff4fc; padding: 15px;">
<span style="font-weight:800;">PRO TIP:</span>
<br><span style="font-style:italic;"> If you are not familiar with Conda yet, please follow the tutorial <a href="https://datascience.101workbook.org/03-SetUpComputingMachine/02C-basic-developer-libraries#conda" target="_blank">Basic Developer Libraries: Conda  ⤴</a> available in the <a href="https://datascience.101workbook.org" target="_blank">Data Science Workbook  ⤴</a></span>.<br><br>
You can <b>install conda</b> following the instructions for regular installation on your operating system: <br>
<b><a href="https://conda.io/projects/conda/en/latest/user-guide/install/windows.html" target="_blank">Windows  ⤴</a></b><br>
<b><a href="https://conda.io/projects/conda/en/latest/user-guide/install/linux.html" target="_blank">Linux  ⤴</a></b><br>
<b><a href="https://conda.io/projects/conda/en/latest/user-guide/install/macos.html" target="_blank">macOS  ⤴</a></b><br>
<i>* If you need a double installation for both ARM and Intel chips on your macOS, please follow the instructions in the tutorial <a href="https://datascience.101workbook.org/03-SetUpComputingMachine/03A-tutorial-installations-on-mac#-install-conda" target="_blank">Installations on MacBook Pro: Install Conda  ⤴</a> available in the <a href="https://datascience.101workbook.org" target="_blank">Data Science Workbook  ⤴</a>.</i><br><br>
Then, you can <b>create a new virtual environment</b> for interactive graphing in the command line: <br>
$ <code>conda create -n graphing python=3.9</code> <br>
activate it:<br>
$ <code>conda activate graphing</code> <br>
and install required libraries:<br>
$ <code>pip install dash</code> <br>
$ <code>pip install dash_bio</code> <br><br>

Once created, the environment can be found using command <code>conda info -e</code> and activated  when needed. Once activated, the additional dependencies can be further installed. To learn more about setting up the environment for interactive graphing, follow the tutorial <a href="https://datascience.101workbook.org/08-DataVisualization/02-GRAPHS/02-PYTHON/01-interactive-graphing-with-python" target="_blank">Interactive Plotting with Python  ⤴</a> available in the <a href="https://datascience.101workbook.org" target="_blank">Data Science Workbook  ⤴</a>.
</div><br>


**Getting started**

First, activate your <a href="https://conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank">Conda  ⤴</a> environment dedicated to <a href="https://datascience.101workbook.org/08-DataVisualization/02-GRAPHS/02-PYTHON/01-interactive-graphing-with-python" target="_blank">Interactive Plotting with Python  ⤴</a>:

```
conda activate graphing
```

![Conda activate](../assets/images/ideogram_conda_activate.png)

Then navigate the file system to the `ideogram` directory in your `{user}.github.io` repository.

```
cd path_to/{user}.github.io/ideogram/
```

![](../assets/images/ideogram_jupyter.png)


Finally, launch `jupyter lab` from the command line and go to the web browser searching for the URL: <a href="http://localhost:8889/lab" target="_blank">http://localhost:8889/lab  ⤴</a>, unless it opens automatically.

```
jupyter lab
```

![Launch Ideogram](../assets/images/ideogram_jupyter.gif)

By default, Jupyter Lab shows on the browser the local file system from the level of the current directory (*where it was started in the terminal*). If the `ideogram-dash.ipynb` notebook does not start automatically in the right-hand panel, click on it twice (*in the browser on the left-hand side*). Then, from the top menu bar, select `Kernel` followed by `Restart kernel and run all cells...`. Scroll down the script and wait until the calculations are complete. Once the output cell appears on the screen, click on the link *http://127.0.0.1:8050/*. In the new tab in your browser, it will open the interactive interface of the Ideogram application. Switch the web tab and adjust options for your inputs and diagram dimensions. Hovering on the options' labels displays the extended instructions and descriptions.

![Ideogram options](../assets/images/ideogram_options.png)

**Use custom bands and annotations**

The `ideogram-dash` application provided in the <a href="https:https://github.com/ISUgenomics" target="_blank">ISUgenomics</a>/<a href="https:https://github.com/ISUgenomics/ideogram_db" target="_blank">ideogram_db  ⤴</a> repository allows you to display a diagram for data stored in any online database. By default, it uses the original **UNPKG** database, available at <a href="https://unpkg.com/ideogram/dist/data/" target="_blank">https://unpkg.com/ideogram/dist/data/  ⤴</a>. There are **bands** and **annotations** subdirectories, respectively. Instead, you can provide a URL to your custom inputs stored in the database hosted on your `{user}.github.io`, created previously in step <a href="https://bioinformaticsworkbook.org/dataVisualization/Plotly/01-ideogram-chromosome-bands.html#upload-data-into-the-database" target="_blank">Upload data into the database  ⤴</a> of this tutorial.

As example, I've created the `aedawid.github.io` repository, available at <a href="https://github.com/aedawid/aedawid.github.io" target="_blank">https://github.com/aedawid/aedawid.github.io  ⤴</a>, which is hosted online via GitHub Pages at <a href="https://aedawid.github.io/" target="_blank">https://aedawid.github.io/  ⤴</a>. The page's index allows entering the `ideogram/database` and displaying the example inputs for bands and annotations.

Let's try to use this data in the interactive `ideogram-dash` web application:

For the option ***Provide URL to online bands*** use:
```
https://aedawid.github.io/ideogram/database/bands/
```

![Select input bands](../assets/images/ideogram_select_bands.png)

*Note that the list of available inputs is updated when a new database is indicated. And when you select a file from that drop-down menu, the list of available chromosomes is updated too.*

For the option ***Provide URL to online annotations*** use:

```
https://aedawid.github.io/ideogram/database/annotations/
```

![Select input annotations](../assets/images/ideogram_select_annots.png)


**Now it's your turn to try your own data!**


### ***B. Use plotly variant***

<span style="color: #ff3870;font-weight: 500;">Section under development...</span>

<!--
### *Convert data structure*

### *Open ideogram in JupyterLab*

-->
