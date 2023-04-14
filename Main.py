import ADN_library as al
import streamlit as st
import Visual_Streamlit as vs
import altair as alt
import pandas as pd



header, ADN_list = al.read_file("Data/sequence.fasta")
df_ADN = al.fromFileToDataframe("Data/sequence.fasta")
vs.head()
vs.bar_chart(df_ADN)

_, df_pourcent_ADN = al.pourcentNucleo(ADN_list)
bar_chart = alt.Chart(df_pourcent_ADN).mark_bar().encode(
    y='Pourcent:Q',
    x='Nucleotide:O',
)

st.altair_chart(bar_chart, use_container_width=True)



