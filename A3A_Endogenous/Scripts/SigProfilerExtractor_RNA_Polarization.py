from SigProfilerExtractor import sigpro as sig
def main_function():    
   # to get input from table format (mutation catalog matrix)
   path_to_example_table = sig.importdata("matrix")
   data = "/project/xiaojiac_1206/SigProfilerExtractor/a3a-polarization-matrix/polarization_matrix_updated-transposed.tsv"
   sig.sigProfilerExtractor("matrix", "a3a-polarization-extractor", data, opportunity_genome="GRCh37", minimum_signatures=1, maximum_signatures=3)
if __name__=="__main__":
   main_function()
