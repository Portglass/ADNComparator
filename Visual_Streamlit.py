import streamlit as st

def head():
    st.title("ADN Comparator")
    activity = ['Intro', 'DNA', 'DotPlot', "About"]
    choice = st.sidebar.selectbox("Select Activity", activity)
    return choice

def bar_chart(df):
    st.bar_chart(df)


