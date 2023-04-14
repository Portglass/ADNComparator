import streamlit as st

def head():
    st.header("ADN Comparator")
    st.subheader("Comparer vos ADN et tirez en maximun d'information")

def bar_chart(df):
    st.bar_chart(df)


