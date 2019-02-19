if __name__ == "__main__":

    #run the test data
    mutFile = "Pipeline/Data/testMutations_p2"
    proteinFolder = "/scratch_ssd/Data/proteins_p2/"

    #prepare the file to submit..
    fasta_file_name = "p2_test_sequences.fasta"
    mutation_sorted_name = "p2_test_mutations_sorted"
    f_file = open(fasta_file_name, 'w')
    m_file = open(mutation_sorted_name, 'w')

    p2_except = ['115298678', '14550407', '38016947', '4502495', '4557729', '45580688', '62414289', '6996021', 'NP_056243.1']

    with open(mutFile,'r') as file:
        currentSeq = None
        i = 0
        for line in file:
                
            line = line.strip()
            if line[0] == ">":
                seqName = line[1:]

                if seqName in p2_except:
                    continue

                #get the fasta sequence for this identifier
                fastaFile = proteinFolder + seqName + '.fasta'
                with open(fastaFile,'r') as file:
                    file.readline()
                    fullSequence = file.read()
                    fullSequence = fullSequence.strip()
                    fullSequence = ''.join(fullSequence.split())


                f_file.write(">%s\n%s\n" % (seqName,fullSequence))

                #write the seqname in mutation file
                m_file.write(">%s\n" % (seqName))

            else:

                #write variant itself to mutation file
                if seqName in p2_except:
                    continue

                #we should also do a check to see if mutation and fasta match on position.
                mut_pos = int(line[1:-1])
                ref_aa = line[0]
                if fullSequence[mut_pos-1] != ref_aa:
                    continue


                m_file.write(line + '\n')





    f_file.close()
    m_file.close()
                