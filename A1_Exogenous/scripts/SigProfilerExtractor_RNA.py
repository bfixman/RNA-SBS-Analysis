from SigProfilerExtractor import sigpro as sig
def main_function():    
   # to get input from table format (mutation catalog matrix)
   data = "/project/xiaojiac_1206/SigProfilerExtractor/a1_matrix/matrix_transposed.tsv"
   sig.sigProfilerExtractor("matrix", "A1_extractor", data, opportunity_genome="GRCh38", minimum_signatures=1, maximum_signatures=3)
if __name__=="__main__":
   main_function()
