from SigProfilerAssignment import Analyzer as Analyze
Analyze.decompose_fit(samples="/project/xiaojiac_1206/SigProfilerExtractor/A3G_exogenous_extractor/CH192/Samples.txt", 
                      output="a3g_exogenous_RNA-SBS",
                      input_type="matrix",
                      signatures="/project/xiaojiac_1206/SigProfilerExtractor/A3G_exogenous_extractor/CH192/Suggested_Solution/CH192_De-Novo_Solution/Signatures/CH192_De-Novo_Signatures.txt",
                      signature_database="/project/xiaojiac_1206/SigProfilerAssignment/Reference/COSMIC_v3.4_RNA-SBS_GRCh37.txt",
                      genome_build="GRCh38")
