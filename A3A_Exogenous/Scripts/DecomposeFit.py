from SigProfilerAssignment import Analyzer as Analyze
Analyze.decompose_fit(samples="/project/xiaojiac_1206/SigProfilerExtractor/A3A_RNA_output/CH192/Samples.txt", 
                      output="A3A_RNA_SBS",
                      input_type="matrix",
                      signatures="/project/xiaojiac_1206/SigProfilerExtractor/A3A_RNA_output/CH192/Suggested_Solution/CH192_De-Novo_Solution/Signatures/CH192_De-Novo_Signatures.txt",
                      signature_database="/project/xiaojiac_1206/SigProfilerAssignment/Reference/COSMIC_v3.4_RNA-SBS_GRCh37.txt",
                      genome_build="GRCh38")
