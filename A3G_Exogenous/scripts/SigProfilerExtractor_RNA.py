from SigProfilerExtractor import sigpro as sig
def main_function():    
   # to get input from table format (mutation catalog matrix)
   path_to_example_table = sig.importdata("matrix")
   data = "/project/xiaojiac_1206/SigProfilerExtractor/A3G_exogenous_matrix/matrix_transposed.tsv"
   sig.sigProfilerExtractor("matrix", "A3G_exogenous_extractor", data, opportunity_genome="GRCh38", minimum_signatures=1, maximum_signatures=3)
if __name__=="__main__":
   main_function()
