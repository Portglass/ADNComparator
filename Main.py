# __________________________________________________________________________________________________
# Import functions
import ADN_library as al # Functions to manipulate ADN sequences
import Visual_Streamlit as vs # Fucntions to design the streamlit page
# __________________________________________________________________________________________________
# Import libraries
import streamlit as st
import matplotlib
matplotlib.use("Agg")
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from io import StringIO
import altair as alt
import pandas as pd
from Bio.pairwise2 import format_alignment

# __________________________________________________________________________________________________
def main():
    choice = vs.head()  # Set up a page selector
    # __________________________________________________________________________________________________
    # Intro Page
    if choice == 'Introduction':
        st.write(":blue[Julie LANGLOIS, Balthazar RAFFY, Martin REVILLON, Line SAIDAN]")
        st.subheader("Vous pouvez réaliser différentes opérations d'analyse et de transcription protéiques sur des "
                 "séquences"
                 " nucléotidiques, à partir de fichiers .fasta ou générées aléatoirement. Vous pouvez également réaliser"
                 " un"
                 " alignement de deux séquences à partir de deux fichiers .fasta.")
    # __________________________________________________________________________________________________
    # Page for .fasta file analysis
    elif choice == "Analyse d'un fichier .fasta":
        st.subheader("DNA details")
        seq_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa"], accept_multiple_files=True)  # Drag&Drop the fasta file into seq_file
        for file in seq_file:  # Use every fasta file uploaded with a for loop
            st.subheader(file.name)  # Print the name of the file
            stringio = StringIO(file.getvalue().decode("utf-8"))  # Decode a file into a string
            sequence = ""  # The string that will receive the DNA sequence
            for record in SeqIO.parse(stringio, 'fasta'):  # get the description and the sequence of a fasta file format
                seq_name = record.description  # Get the description from the fasta format
                sequence = str(record.seq)  # Get the sequence as a string
            st.subheader(seq_name)  # print the sequence description
            st.write(al.summaryADNText(list(sequence)))  # Print the numb of nucleotides and if the sequence is correct
            # ___________________________________________________________________________________________________
            df_dna = al.fromStringToDataframe(list(sequence))  # convert the sequence as a df to display it as barchart
            st.subheader(':blue[DNA Sequence]')
            vs.bar_chart(df_dna)  # Display the DNA sequence as a barchart
            # ___________________________________________________________________________________________________
            st.subheader(':blue[Frequency of each base]')
            _, df_pourcent_ADN = al.pourcentNucleo(list(sequence))  # Get the frequency of each base as a dataframe
            bar_chart = alt.Chart(df_pourcent_ADN).mark_bar().encode(
                y='Pourcent:Q',
                x='Nucleotide:O',
            )
            st.altair_chart(bar_chart, use_container_width=True)
            # Display the barchart of each frequency
            # __________________________________________________________________________________________________
            comp_seq = al.complementaire(sequence)  # Get the complementary sequence
            rev_seq = al.transposer(comp_seq)  # Reverse the complementary sequence
            df_rev_dna = al.fromStringToDataframe(list(rev_seq))  # Transform rev_seq into a dataframe
            st.subheader(':blue[ Complementary DNA Sequence]')
            vs.bar_chart(df_rev_dna)  # Display the DNA sequence as a barchart
            # __________________________________________________________________________________________________
            st.write("________________________")
            dic = al.seq_dic(sequence)  # Create dictionary with every 4 nucleotides patterns
            df_dic = pd.DataFrame.from_dict(dic,orient='index')  # Transform the dictionary into dataframe
            st.subheader(':blue[ Motifs nucléotidiques]')
            st.write('Il y a ' + str(len(dic)) + ' motifs distincts de 4 nucléotides dans la séquence.')
            vs.bar_chart(df_dic)  # Display the dataframe
            # __________________________________________________________________________________________________
            st.subheader(":blue[Chercher la fréquence d'un motif nucléotidiques]")
            pattern = st.text_input('Veuillez entrer un motif:')  # User input a pattern
            if pattern == '':  # Display nothing when no pattern is entered
                st.write("")
            elif not al.isCorrectADN(pattern):  # Check if the pattern is a valid sequence
                st.write(':red[Le motif entré est incorrect !]')
            elif len(pattern) > len(sequence):  # Check that the pattern is not longer than the sequence
                st.write(':red[Le motif entré est trop long !]')
            else:
                freq, position = al.freq_pattern(sequence, pattern)  # Count frequency of a correct pattern
                if freq != 0:
                    st.write('Il y a ' + str(freq) + ' occurences de ce motif aux positions suivantes:' + str(position))
                else:
                    st.write("Ce motif n'est pas présent dans la séquence !")
            # __________________________________________________________________________________________________
            st.write("________________________")
            st.subheader(":blue[Codons STOP]")
            f_TAA, p_TAA = al.freq_pattern(sequence,'TAA')
            f_TAG, p_TAG = al.freq_pattern(sequence, 'TAG')
            f_TGA, p_TGA = al.freq_pattern(sequence, 'TGA')
            st.write('Il y a ' + str(f_TGA + f_TAA + f_TAG) + ' codons STOP aux positions suivantes:' + str(p_TGA + p_TAG + p_TAA))
            st.write("________________________")
            # __________________________________________________________________________________________________
            st.subheader(":blue[Traduction en protéine]")
            orf1 = Seq(sequence)
            orf2 = Seq(sequence[1:])
            orf3 = Seq(sequence[2:])
            prot1 = orf1.translate()
            prot2 = orf2.translate()
            prot3 = orf3.translate()
            # ________________________
            st.write("\n:blue[1er cadre de lecture]")
            st.write(str(prot1))
            st.write(':blue[Frequency of each amino acid (%)]')
            dic_prot1 = ProteinAnalysis(str(prot1)).get_amino_acids_percent()
            df_dic_prot1 = pd.DataFrame.from_dict(dic_prot1, orient='index')  # Transform the dictionary into dataframe
            vs.bar_chart(df_dic_prot1.mul(100))  # Display the dataframe
            st.write(":blue[Chercher la fréquence d'un motif d'acides aminés]")
            pattern_prot1= st.text_input('Veuillez entrer un motif pour le 1er cadre de lecture:')  # User input a pattern
            if pattern_prot1 == '':  # Display nothing when no pattern is entered
                st.write("")
            elif not al.isCorrectaa(pattern_prot1):  # Check if the pattern is a valid sequence
                st.write(':red[Le motif entré est incorrect !]')
            elif len(pattern_prot1) > len(prot1):  # Check that the pattern is not longer than the sequence
                st.write(':red[Le motif entré est trop long !]')
            else:
                freq_1, position_1 = al.freq_pattern(prot1, pattern_prot1)  # Count frequency of a correct pattern
                if freq_1 != 0:
                    st.write('Il y a ' + str(freq_1) + ' occurences de ce motif aux positions suivantes:' + str(position_1))
                else:
                    st.write("Ce motif n'est pas présent dans la séquence !")
            # ________________________
            st.write("________________________")
            st.write("\n:blue[2ème cadre de lecture]")
            st.write(str(prot2))
            st.write(':blue[Frequency of each amino acid (%)]')
            dic_prot2 = ProteinAnalysis(str(prot2)).get_amino_acids_percent()
            df_dic_prot2 = pd.DataFrame.from_dict(dic_prot2, orient='index')  # Transform the dictionary into dataframe
            vs.bar_chart(df_dic_prot2.mul(100))  # Display the dataframe
            st.write(":blue[Chercher la fréquence d'un motif d'acides aminés]")
            pattern_prot2 = st.text_input(
                'Veuillez entrer un motif pour le 2ème cadre de lecture:')  # User input a pattern
            if pattern_prot2 == '':  # Display nothing when no pattern is entered
                st.write("")
            elif not al.isCorrectaa(pattern_prot2):  # Check if the pattern is a valid sequence
                st.write(':red[Le motif entré est incorrect !]')
            elif len(pattern_prot2) > len(prot2):  # Check that the pattern is not longer than the sequence
                st.write(':red[Le motif entré est trop long !]')
            else:
                freq_2, position_2 = al.freq_pattern(prot2, pattern_prot2)  # Count frequency of a correct pattern
                if freq_2 != 0:
                    st.write(
                        'Il y a ' + str(freq_2) + ' occurences de ce motif aux positions suivantes:' + str(position_2))
                else:
                    st.write("Ce motif n'est pas présent dans la séquence !")
            # _________________________
            st.write("________________________")
            st.write("\n:blue[3ème cadre de lecture]")
            st.write(str(prot3))
            st.write(':blue[Frequency of each amino acid (%)]')
            dic_prot3 = ProteinAnalysis(str(prot3)).get_amino_acids_percent()
            df_dic_prot3 = pd.DataFrame.from_dict(dic_prot3, orient='index')  # Transform the dictionary into dataframe
            vs.bar_chart(df_dic_prot3.mul(100))  # Display the dataframe
            st.write(":blue[Chercher la fréquence d'un motif d'acides aminés]")
            pattern_prot3 = st.text_input('Veuillez entrer un motif pour le 3ème cadre de lecture:')  # User input a pattern
            if pattern_prot3 == '':  # Display nothing when no pattern is entered
                st.write("")
            elif not al.isCorrectaa(pattern_prot3):  # Check if the pattern is a valid sequence
                st.write(':red[Le motif entré est incorrect !]')
            elif len(pattern_prot3) > len(prot3):  # Check that the pattern is not longer than the sequence
                st.write(':red[Le motif entré est trop long !]')
            else:
                freq_3, position_3 = al.freq_pattern(prot3, pattern_prot3)  # Count frequency of a correct pattern
                if freq_3 != 0:
                    st.write(
                        'Il y a ' + str(freq_3) + ' occurences de ce motif aux positions suivantes:' + str(position_3))
                else:
                    st.write("Ce motif n'est pas présent dans la séquence !")
    # __________________________________________________________________________________________________
    # Page to analyse a random nucleotides sequence
    elif choice == "Analyse d'une séquence aléatoire":
        st.subheader("Analyse d'une séquence aléatoire")
        # __________________________________________________________________________________________________
        st.subheader(":blue[Générer un séquence aléatoire]")
        size_seq = st.number_input('Taille de la séquence aléatoire:', step=1)
        if size_seq != 0:
            sequence = al.random_seq(int(size_seq))
            df_random_seq = al.fromStringToDataframe(list(sequence))  # Transform rev_seq into a dataframe
            st.subheader(':blue[Séquence Aléatoire]')
            vs.bar_chart(df_random_seq)
            # ___________________________________________________________________________________________________
            st.subheader(':blue[Frequency of each base]')
            _, df_pourcent_dna = al.pourcentNucleo(list(sequence))  # Get the frequency of each base as a dataframe
            bar_chart = alt.Chart(df_pourcent_dna).mark_bar().encode(
                y='Pourcent:Q',
                x='Nucleotide:O',
            )
            st.altair_chart(bar_chart, use_container_width=True)  # Display the barchart of each frequency
            # __________________________________________________________________________________________________
            comp_seq = al.complementaire(sequence)  # Get the complementary sequence
            rev_seq = al.transposer(comp_seq)  # Reverse the complementary sequence
            df_rev_dna = al.fromStringToDataframe(list(rev_seq))  # Transform rev_seq into a dataframe
            st.subheader(':blue[ Complementary DNA Sequence]')
            vs.bar_chart(df_rev_dna)  # Display the DNA sequence as a barchart
            # __________________________________________________________________________________________________
            dic = al.seq_dic(sequence)  # Create dictionary with every 4 nucleotides patterns
            df_dic = pd.DataFrame.from_dict(dic, orient='index')  # Transform the dictionary into dataframe
            st.subheader(':blue[ Motifs nucléotidiques]')
            st.write('Il y a ' + str(len(dic)) + ' motifs distincts de 4 nucléotides dans la séquence.')
            vs.bar_chart(df_dic)  # Display the dataframe
            # __________________________________________________________________________________________________
            st.subheader(":blue[Chercher la fréquence d'un motif nucléotidiques]")
            pattern = st.text_input('Veuillez entrer un motif:')  # User input a pattern
            if pattern == '':  # Display nothing when no pattern is entered
                st.write("")
            elif not al.isCorrectADN(pattern):  # Check if the pattern is a valid sequence
                st.write(':red[Le motif entré est incorrect !]')
            elif len(pattern) > len(sequence):  # Check that the pattern is not longer than the sequence
                st.write(':red[Le motif entré est trop long !]')
            else:
                freq, position = al.freq_pattern(sequence, pattern)  # Count frequency of a correct pattern
                if freq != 0:
                    st.write('Il y a ' + str(freq) + ' occurences de ce motif aux positions suivantes:' + str(position))
                else:
                    st.write("Ce motif n'est pas présent dans la séquence !")
            # __________________________________________________________________________________________________
            st.subheader(":blue[Codons STOP]")
            f_TAA, p_TAA = al.freq_pattern(sequence, 'TAA')
            f_TAG, p_TAG = al.freq_pattern(sequence, 'TAG')
            f_TGA, p_TGA = al.freq_pattern(sequence, 'TGA')
            st.write('Il y a ' + str(f_TGA + f_TAA + f_TAG) + ' codons STOP aux positions suivantes:' + str(
                p_TGA + p_TAG + p_TAA))
            # __________________________________________________________________________________________________
            st.subheader(":blue[Traduction en protéine]")
            orf1 = Seq(sequence)
            orf2 = Seq(sequence[1:])
            orf3 = Seq(sequence[2:])
            prot1 = orf1.translate()
            prot2 = orf2.translate()
            prot3 = orf3.translate()
            st.write(":blue[1er cadre de lecture]")
            st.write(str(prot1))
            st.subheader(':blue[Frequency of each amino acid]')
            _, df_pourcent_prot = al.pourcentNucleo(list(prot1))  # Get the frequency of each base as a dataframe
            bar_chart = alt.Chart(df_pourcent_prot).mark_bar().encode(
                y='Pourcent:Q',
                x='Nucleotide:O',
            )
            st.altair_chart(bar_chart, use_container_width=True)
            st.write(":blue[2ème cadre de lecture]")
            st.write(str(prot2))
            st.write(":blue[3ème cadre de lecture]")
            st.write(str(prot3))
            st.write("\n:blue[1er cadre de lecture]")
            st.write(str(prot1))
            st.write(':blue[Frequency of each amino acid (%)]')
            dic_prot1 = ProteinAnalysis(str(prot1)).get_amino_acids_percent()
            df_dic_prot1 = pd.DataFrame.from_dict(dic_prot1, orient='index')  # Transform the dictionary into dataframe
            vs.bar_chart(df_dic_prot1.mul(100))  # Display the dataframe
            st.write(":blue[Chercher la fréquence d'un motif d'acides aminés]")
            pattern_prot1 = st.text_input(
                'Veuillez entrer un motif pour le 1er cadre de lecture:')  # User input a pattern
            if pattern_prot1 == '':  # Display nothing when no pattern is entered
                st.write("")
            elif not al.isCorrectaa(pattern_prot1):  # Check if the pattern is a valid sequence
                st.write(':red[Le motif entré est incorrect !]')
            elif len(pattern_prot1) > len(prot1):  # Check that the pattern is not longer than the sequence
                st.write(':red[Le motif entré est trop long !]')
            else:
                freq_1, position_1 = al.freq_pattern(prot1, pattern_prot1)  # Count frequency of a correct pattern
                if freq_1 != 0:
                    st.write(
                        'Il y a ' + str(freq_1) + ' occurences de ce motif aux positions suivantes:' + str(position_1))
                else:
                    st.write("Ce motif n'est pas présent dans la séquence !")
            # ________________________
            st.write("________________________")
            st.write("\n:blue[2ème cadre de lecture]")
            st.write(str(prot2))
            st.write(':blue[Frequency of each amino acid (%)]')
            dic_prot2 = ProteinAnalysis(str(prot2)).get_amino_acids_percent()
            df_dic_prot2 = pd.DataFrame.from_dict(dic_prot2, orient='index')  # Transform the dictionary into dataframe
            vs.bar_chart(df_dic_prot2.mul(100))  # Display the dataframe
            st.write(":blue[Chercher la fréquence d'un motif d'acides aminés]")
            pattern_prot2 = st.text_input(
                'Veuillez entrer un motif pour le 2ème cadre de lecture:')  # User input a pattern
            if pattern_prot2 == '':  # Display nothing when no pattern is entered
                st.write("")
            elif not al.isCorrectaa(pattern_prot2):  # Check if the pattern is a valid sequence
                st.write(':red[Le motif entré est incorrect !]')
            elif len(pattern_prot2) > len(prot2):  # Check that the pattern is not longer than the sequence
                st.write(':red[Le motif entré est trop long !]')
            else:
                freq_2, position_2 = al.freq_pattern(prot2, pattern_prot2)  # Count frequency of a correct pattern
                if freq_2 != 0:
                    st.write(
                        'Il y a ' + str(freq_2) + ' occurences de ce motif aux positions suivantes:' + str(position_2))
                else:
                    st.write("Ce motif n'est pas présent dans la séquence !")
            # _________________________
            st.write("________________________")
            st.write("\n:blue[3ème cadre de lecture]")
            st.write(str(prot3))
            st.write(':blue[Frequency of each amino acid (%)]')
            dic_prot3 = ProteinAnalysis(str(prot3)).get_amino_acids_percent()
            df_dic_prot3 = pd.DataFrame.from_dict(dic_prot3, orient='index')  # Transform the dictionary into dataframe
            vs.bar_chart(df_dic_prot3.mul(100))  # Display the dataframe
            st.write(":blue[Chercher la fréquence d'un motif d'acides aminés]")
            pattern_prot3 = st.text_input(
                'Veuillez entrer un motif pour le 3ème cadre de lecture:')  # User input a pattern
            if pattern_prot3 == '':  # Display nothing when no pattern is entered
                st.write("")
            elif not al.isCorrectaa(pattern_prot3):  # Check if the pattern is a valid sequence
                st.write(':red[Le motif entré est incorrect !]')
            elif len(pattern_prot3) > len(prot3):  # Check that the pattern is not longer than the sequence
                st.write(':red[Le motif entré est trop long !]')
            else:
                freq_3, position_3 = al.freq_pattern(prot3, pattern_prot3)  # Count frequency of a correct pattern
                if freq_3 != 0:
                    st.write(
                        'Il y a ' + str(freq_3) + ' occurences de ce motif aux positions suivantes:' + str(position_3))
                else:
                    st.write("Ce motif n'est pas présent dans la séquence !")

    # __________________________________________________________________________________________________
    # Align two .fasta file
    elif choice == 'Alignement de séquences':
        st.subheader("Comparaison entre deux ADN")
        st.subheader(":blue[Entrez la 1ère séquence]")
        seq_file1 = st.file_uploader("Upload first FASTA File", type=["fasta", "fa"],
                                    accept_multiple_files=True)  # Drag&Drop the fasta file into seq_file
        for file in seq_file1:  # Use every fasta file uploaded with a for loop
            st.write(file.name)  # Print the name of the file
            stringio = StringIO(file.getvalue().decode("utf-8"))  # Decode a file into a string
            sequence1 = ""  # The string that will receive the DNA sequence
            for record in SeqIO.parse(stringio, 'fasta'):  # get the description and the sequence of a fasta file format
                seq_name1 = record.description  # Get the description from the fasta format
                sequence1 = str(record.seq)  # Get the sequence as a string
            st.subheader(seq_name1)  # print the sequence description
            st.write(al.summaryADNText(list(sequence1)))  # Print the numb of nucleotides and if the sequence is correct
            # list_seq.append(list(sequence))  # Add the current sequence to the list
            # ___________________________________________________________________________________________________
            df_dna1 = al.fromStringToDataframe(list(sequence1))  # convert the sequence as a df to display it as barchart
            st.subheader(':blue[DNA Sequence n°1]')
            vs.bar_chart(df_dna1)  # Display the DNA sequence as a barchart
            # ___________________________________________________________________________________________________
        st.write("___________________________________")
        st.subheader(":blue[Entrez la 2ème séquence]")
        seq_file2 = st.file_uploader("Upload second FASTA File", type=["fasta", "fa"],
                                     accept_multiple_files=True)  # Drag&Drop the fasta file into seq_file
        for file in seq_file2:  # Use every fasta file uploaded with a for loop
            st.write(file.name)  # Print the name of the file
            stringio = StringIO(file.getvalue().decode("utf-8"))  # Decode a file into a string
            sequence2 = ""  # The string that will receive the DNA sequence
            for record in SeqIO.parse(stringio,
                                      'fasta'):  # get the description and the sequence of a fasta file format
                seq_name2 = record.description  # Get the description from the fasta format
                sequence2 = str(record.seq)  # Get the sequence as a string
            st.subheader(seq_name2)  # print the sequence description
            st.write(
                al.summaryADNText(list(sequence2)))  # Print the numb of nucleotides and if the sequence is correct
            # list_seq.append(list(sequence))  # Add the current sequence to the list
            # ___________________________________________________________________________________________________
            df_dna2 = al.fromStringToDataframe(
                list(sequence2))  # convert the sequence as a df to display it as barchart
            st.subheader(':blue[DNA Sequence n°1]')
            vs.bar_chart(df_dna2)  # Display the DNA sequence as a barchart
        # __________________________________________________________________________________________________
            st.write("___________________________________")
            if st.button ("Réaliser l'alignement"):
                if sequence2 == sequence1 :
                    st.write(":red[Les deux séquences sont identiques!]")
                else:
                    align = al.align_seq(sequence1,sequence2)
                    for a in align:
                        st.write(format_alignment(*a))



    # __________________________________________________________________________________________________


if __name__ == '__main__':
    main()
