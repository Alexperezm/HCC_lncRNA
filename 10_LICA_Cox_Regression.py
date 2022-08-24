# -*- coding: utf-8 -*-
"""10_LICA_COX_Regression.ipynb

10.Cox Regression - LICA

Author: Alexperezm | Master's End of Degree Project - 2021-2022

Objective: Study and development of predictive models of survival based on several clinial variables (Gender, age, alcohol intake, Hepatitis C, Hepatitis B, Fibrosis stage and Edmonson Grade).

Loading the libraries:
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

!pip install lifelines
#import lifelines

from lifelines import KaplanMeierFitter


df = pd.read_excel("/content/drive/MyDrive/TFM/Metadata_LICA.xlsx")
df.head()

"""The dataset must be processed in order to adeccuate the variables for the modelling."""

Fact_var = df[["Alcohol intake", "Hepatitis B","Hepatitis C", "Metabolic syndrome","Vascular invasion"]]

Fact_var

cat = pd.Categorical(Fact_var["Alcohol intake"],categories=["no","yes"])
codes, uniques = pd.factorize(cat)
Fact_var["Alcohol intake"] = codes

cat = pd.Categorical(Fact_var["Hepatitis B"],categories=["no","yes"])
codes, uniques = pd.factorize(cat)
Fact_var["Hepatitis B"] = codes

cat = pd.Categorical(Fact_var["Hepatitis C"],categories=["no","yes"])
codes, uniques = pd.factorize(cat)
Fact_var["Hepatitis C"] = codes

cat = pd.Categorical(Fact_var["Metabolic syndrome"],categories=["no","yes"])
codes, uniques = pd.factorize(cat)
Fact_var["Metabolic syndrome"] = codes

cat = pd.Categorical(Fact_var["Vascular invasion"],categories=["no","yes"])
codes, uniques = pd.factorize(cat)
Fact_var["Vascular invasion"] = codes

Fact_var

cat = pd.Categorical(df["Gender"], categories= ["M", "F"])
codes, uniques = pd.factorize(cat)
Fact_var["Gender"] = codes
print(uniques)
#Fact_var['Gender']= df["Gender"] == "M"
#Fact_var['Gender']= Fact_var['Gender']*1

cat = pd.Categorical(df["Fibrosis stage"], categories= ["F0-F1", "F2-F3","F4"])
codes, uniques = pd.factorize(cat)
Fact_var["Fibrosis stage"] = codes
print(uniques)

cat = pd.Categorical(df["Edmonson grade"], categories= ["I-II", "III-IV"])
codes, uniques = pd.factorize(cat)
Fact_var["Edmonson grade"] = codes
print(uniques)

#Fact_var['Edmonson grade']= df["Edmonson grade"] == "III-IV"
#Fact_var['Edmonson grade']= Fact_var['Edmonson grade']*1

Fact_var["Age"] = df["Age"]
Fact_var["Largest nodule diameter (mm)"] = df["Largest nodule diameter (mm)"]
Fact_var["Last survival (days)"] = df ["Last survival (days)"]
Fact_var["Vital status"] = df["Vital status"]

Fact_var

"""Review: all variables are codified correctly."""

Fact_var.dtypes

Fact_var.isnull().sum()

Fact_var = Fact_var.dropna()

Fact_var.isnull().sum()

T = Fact_var["Last survival (days)"]
E = Fact_var["Vital status"]
plt.hist(T, bins = 50)
plt.show()

"""# Dummy Variables:

For those variables displaying several stages are made dummy variables, in order to simplify the modeling.
"""

dummies_ecog = pd.get_dummies(Fact_var["Fibrosis stage"], prefix = 'Fibrosis stage')
dummies_ecog.head(4)

dummies_ecog = dummies_ecog[["Fibrosis stage_0", "Fibrosis stage_1", "Fibrosis stage_2"]]
Fact_var = pd.concat([Fact_var, dummies_ecog], axis = 1)
Fact_var = Fact_var.drop("Fibrosis stage", axis = 1)
Fact_var.head()

dummies_ecog = pd.get_dummies(Fact_var["Edmonson grade"], prefix = 'Edmonson grade')
dummies_ecog.head(4)

dummies_ecog = dummies_ecog[["Edmonson grade_-1", "Edmonson grade_0", "Edmonson grade_1"]]
Fact_var = pd.concat([Fact_var, dummies_ecog], axis = 1)
Fact_var = Fact_var.drop("Edmonson grade", axis = 1)
Fact_var.head()

"""# Non parametric testing:

Development of the survival plots:
"""

kmf = KaplanMeierFitter()
kmf.fit(durations = T, event_observed = E)
kmf.plot_survival_function()

kmf.survival_function_.plot()
plt.title('Survival function')

kmf.plot_cumulative_density()

from lifelines.utils import median_survival_times
median_ = kmf.median_survival_time_
median_confidence_interval_ = median_survival_times(kmf.confidence_interval_)
print(median_)
print(median_confidence_interval_)

"""# Cox Proportional Hazard Model (Semi-Parametric)

Correlation matrix:

Based on Cox Proportional Hazard Model, is done the correlation matrix and afterwards the Kaplan Meier for each variable. The whole procedure is thoroughly explained in the memoir.
"""

from lifelines import CoxPHFitter
cph = CoxPHFitter()
FIT=CoxPHFitter(penalizer=0.1).fit(Fact_var, duration_col = 'Last survival (days)', event_col = 'Vital status')
FIT.print_summary()

FIT.plot()

matrix = Fact_var.corr()
print(matrix)

ax = plt.subplot(111)
m = (Fact_var["Gender"] == 0)

kmf_male = KaplanMeierFitter()
ax = kmf_male.fit(durations = T[m], event_observed = E[m], label = "Hombre").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)
ax.set_xlabel('Días', size = 10)
ax.set_ylabel('Ratio Supervivencia', size = 10)
np.arange(0, 1, step=0.2)
ax.set_xticks([0,365,730,1095,1460,1825]) #np.arange(0, 1825, step=365)

kmf_female = KaplanMeierFitter()
ax = kmf_female.fit(durations = T[~m], event_observed = E[~m], label = "Mujer").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)
#time_limit = 1825

from lifelines.plotting import add_at_risk_counts
add_at_risk_counts(kmf_female, kmf_male, ax=ax,xticks=[365,730,1095,1460,1825])


#plt.
#plt.tight_layout()
plt.xlabel("Supervivencia - días")
plt.title("Supervivencia según género")
#plt.text(6, 9.5, 'Here we go')

from lifelines.statistics import survival_difference_at_fixed_point_in_time_test
point_in_time = 1825
results = survival_difference_at_fixed_point_in_time_test(point_in_time, kmf_female, kmf_male)
results.print_summary()

ax = plt.subplot(111)
m = (Fact_var["Alcohol intake"] == 0)

kmf_abst = KaplanMeierFitter()
ax = kmf_abst.fit(durations = T[m], event_observed = E[m], label = "Abstinencia").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)
ax.set_xlabel('Días', size = 10)
ax.set_ylabel('Ratio Supervivencia', size = 10)
np.arange(0, 1, step=0.2)
ax.set_xticks([0,365,730,1095,1460,1825])

kmf_cons = KaplanMeierFitter()
ax = kmf_cons.fit(T[~m], event_observed = E[~m], label = "Consumo").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)

from lifelines.plotting import add_at_risk_counts
add_at_risk_counts(kmf_cons, kmf_abst, ax=ax,xticks=[365,730,1095,1460,1825])

#plt.tight_layout()
plt.xlabel("Supervivencia - días")
plt.title("Supervivencia según consumo de alcohol")

from lifelines.statistics import survival_difference_at_fixed_point_in_time_test
point_in_time = 1825
results = survival_difference_at_fixed_point_in_time_test(point_in_time, kmf_abst, kmf_cons)
results.print_summary()

from lifelines.statistics import survival_difference_at_fixed_point_in_time_test
point_in_time = 1460
results = survival_difference_at_fixed_point_in_time_test(point_in_time, kmf_abst, kmf_cons)
results.print_summary()

ax = plt.subplot(111)
m = (Fact_var["Hepatitis C"] == 0)

kmf_healthy = KaplanMeierFitter()
ax = kmf_healthy.fit(durations = T[m], event_observed = E[m], label = "Sano").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)
ax.set_xlabel('Días', size = 10)
ax.set_ylabel('Ratio Supervivencia', size = 10)
np.arange(0, 1, step=0.2)
ax.set_xticks([0,365,730,1095,1460,1825]) 

kmf_hepC = KaplanMeierFitter()
ax = kmf_hepC.fit(T[~m], event_observed = E[~m], label = "Infectado").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)

from lifelines.plotting import add_at_risk_counts
add_at_risk_counts(kmf_hepC, kmf_healthy, ax=ax, xticks=[365,730,1095,1460,1825])

#plt.tight_layout()
plt.xlabel("Supervivencia - días")
plt.title("Supervivencia según Hepatitis C")

from lifelines.statistics import survival_difference_at_fixed_point_in_time_test
point_in_time = 1825
results = survival_difference_at_fixed_point_in_time_test(point_in_time, kmf_healthy, kmf_hepC)
results.print_summary()

ax = plt.subplot(111)
m = (Fact_var["Fibrosis stage_0"] == 1)
kmf_F0_F1 = KaplanMeierFitter()
ax = kmf_F0_F1.fit(durations = T[m], event_observed = E[m], label = "F0-F1").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)
ax.set_xlabel('Días', size = 10)
ax.set_ylabel('Ratio Supervivencia', size = 10)
np.arange(0, 1, step=0.2)
ax.set_xticks([0,365,730,1095,1460,1825]) 

n = (Fact_var["Fibrosis stage_1"] == 1)
kmf_F2_F3 = KaplanMeierFitter()
ax = kmf_F2_F3.fit(T[n], event_observed = E[n], label = "F2-F3").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)

o = (Fact_var["Fibrosis stage_2"] == 1)
kmf_F4 = KaplanMeierFitter()
ax = kmf_F4.fit(T[o], event_observed = E[o], label = "F4").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)

from lifelines.plotting import add_at_risk_counts
add_at_risk_counts(kmf_F4, kmf_F2_F3, ax=ax, xticks=[365,730,1095,1460,1825]) #,kmf_F0_F1

#plt.tight_layout()
plt.xlabel("Supervivencia - días")
plt.title("Supervivencia según grado de fibrosis")

from lifelines.statistics import survival_difference_at_fixed_point_in_time_test
point_in_time = 1825
results = survival_difference_at_fixed_point_in_time_test(point_in_time, kmf_F0_F1, kmf_F2_F3)
results.print_summary()
results = survival_difference_at_fixed_point_in_time_test(point_in_time, kmf_F0_F1, kmf_F4)
results.print_summary()
results = survival_difference_at_fixed_point_in_time_test(point_in_time, kmf_F2_F3,kmf_F4)
results.print_summary()

ax = plt.subplot(111)
m = (Fact_var["Edmonson grade_0"] == 1)
edm_I_II = KaplanMeierFitter()
ax = edm_I_II.fit(durations = T[m], event_observed = E[m], label = "I-II").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)
ax.set_xlabel('Días', size = 10)
ax.set_ylabel('Ratio Supervivencia', size = 10)
np.arange(0, 1, step=0.2)
ax.set_xticks([0,365,730,1095,1460,1825]) 

n = (Fact_var["Edmonson grade_1"] == 1)
edm_III_IV = KaplanMeierFitter()
ax = edm_III_IV.fit(T[n], event_observed = E[n], label = "III-IV").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)


from lifelines.plotting import add_at_risk_counts
add_at_risk_counts(edm_III_IV, edm_I_II, ax=ax, xticks=[365,730,1095,1460,1825])

#plt.tight_layout()
plt.xlabel("Supervivencia - días")
plt.title("Supervivencia según Grado de Edmonson")

from lifelines.statistics import survival_difference_at_fixed_point_in_time_test
point_in_time = 1825
results = survival_difference_at_fixed_point_in_time_test(point_in_time, edm_I_II, edm_III_IV)
results.print_summary()

ax = plt.subplot(111)
m = (Fact_var["Hepatitis B"] == 1)
kmf_healthy = KaplanMeierFitter()
ax = kmf_healthy.fit(durations = T[m], event_observed = E[m], label = "Sano").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)
ax.set_xlabel('Días', size = 10)
ax.set_ylabel('Ratio Supervivencia', size = 10)
np.arange(0, 1, step=0.2)
ax.set_xticks([0,365,730,1095,1460,1825]) 

kmf_hepC = KaplanMeierFitter()
ax = kmf_hepC.fit(T[~m], event_observed = E[~m], label = "Infectado").plot_survival_function(ax=ax)
ax.set_xlim(0,1825)
ax.set_xlabel('Días', size = 10)
ax.set_ylabel('Ratio Supervivencia', size = 10)
np.arange(0, 1, step=0.2)
ax.set_xticks([0,365,730,1095,1460,1825]) 

from lifelines.plotting import add_at_risk_counts
add_at_risk_counts(kmf_hepC, kmf_healthy, ax=ax)

#plt.tight_layout()
plt.xlabel("Supervivencia - días")
plt.title("Supervivencia según Hepatitis B")

from lifelines.statistics import survival_difference_at_fixed_point_in_time_test
point_in_time = 1825
results = survival_difference_at_fixed_point_in_time_test(point_in_time, kmf_healthy, kmf_hepC)
results.print_summary()
