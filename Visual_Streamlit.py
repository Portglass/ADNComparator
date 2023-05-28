# Contains streamlit elements for display
import streamlit as st


def head():
    st.title("Projet de Bioinformatique")
    activity = ['Introduction', "Analyse d'un fichier .fasta", "Analyse d'une séquence aléatoire", 'Alignement de séquences']
    choice = st.sidebar.selectbox("Select Activity", activity)
    return choice


def bar_chart(df):
    st.bar_chart(df, y='')


