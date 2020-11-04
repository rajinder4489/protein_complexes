########################

c("CPX-1027", "CPX-1032", "CPX-104", "CPX-1063", "CPX-1064", "CPX-1066", "CPX-114", "CPX-1147", "CPX-1159", "CPX-1163", "CPX-1164", "CPX-117", "CPX-1194", "CPX-1195", "CPX-1196", "CPX-1199", "CPX-1201", "CPX-1202", "CPX-1203", "CPX-1204", "CPX-1207", "CPX-1209", "CPX-1212", "CPX-1216", "CPX-1219", "CPX-1222", "CPX-1225", "CPX-1226", "CPX-130", "CPX-1313", "CPX-1427", "CPX-148", "CPX-150", "CPX-1724", "CPX-1725", "CPX-1727", "CPX-1728", "CPX-1736", "CPX-175", "CPX-1751", "CPX-1752", "CPX-1770", "CPX-1772", "CPX-1781", "CPX-1782", "CPX-1783", "CPX-1784", "CPX-1796", "CPX-1798", "CPX-1799", "CPX-1801", "CPX-1802", "CPX-1803", "CPX-1804", "CPX-1815", "CPX-1816", "CPX-1817", "CPX-1818", "CPX-1819", "CPX-1820", "CPX-1821", "CPX-1822", "CPX-1826", "CPX-1828", "CPX-1829", "CPX-1842", "CPX-1843", "CPX-1847", "CPX-1910", "CPX-1917", "CPX-1942", "CPX-195", "CPX-1969", "CPX-1971", "CPX-1972", "CPX-2004", "CPX-2005", "CPX-2007", "CPX-2008", "CPX-2011", "CPX-2014", "CPX-2015", "CPX-2016", "CPX-2097", "CPX-2108", "CPX-2159", "CPX-2173", "CPX-2188", "CPX-2201", "CPX-222", "CPX-2375", "CPX-256", "CPX-262", "CPX-2919", "CPX-294", "CPX-2952", "CPX-2953", "CPX-2955", "CPX-3073", "CPX-308", "CPX-312", "CPX-3138", "CPX-3141", "CPX-316", "CPX-3186", "CPX-3199", "CPX-3200", "CPX-3201", "CPX-3227", "CPX-3229", "CPX-3230", "CPX-3232", "CPX-3239", "CPX-3256", "CPX-3263", "CPX-3278", "CPX-329", "CPX-3322", "CPX-354", "CPX-356", "CPX-358", "CPX-3622", "CPX-3624", "CPX-364", "CPX-3946", "CPX-4082", "CPX-4084", "CPX-414", "CPX-4141", "CPX-4142", "CPX-4143", "CPX-4144", "CPX-419", "CPX-4203", "CPX-4206", "CPX-4207", "CPX-4223", "CPX-4224", "CPX-4225", "CPX-4226", "CPX-439", "CPX-442", "CPX-459", "CPX-462", "CPX-467", "CPX-469", "CPX-475", "CPX-483", "CPX-485", "CPX-494", "CPX-497", "CPX-5001", "CPX-502", "CPX-5021", "CPX-5022", "CPX-5029", "CPX-5043", "CPX-5050", "CPX-514", "CPX-5151", "CPX-519", "CPX-52", "CPX-54", "CPX-56", "CPX-5642", "CPX-5643", "CPX-5644", "CPX-5645", "CPX-5646", "CPX-5661", "CPX-5662", "CPX-5663", "CPX-5664", "CPX-578", "CPX-5786", "CPX-5834", "CPX-5840", "CPX-5841", "CPX-5842", "CPX-5843", "CPX-5846", "CPX-614", "CPX-631", "CPX-648", "CPX-654", "CPX-688", "CPX-696", "CPX-709", "CPX-745", "CPX-757", "CPX-77", "CPX-8", "CPX-80", "CPX-880", "CPX-9", "CPX-922", "CPX-955", "CPX-974", "CPX-978", "CPX-979", "CPX-98", "CPX-985", "CPX-99")

test <- complex_file[complex_file$Complex.ac %in% sel_cids, "Identifiers..and.stoichiometry..of.molecules.in.complex"]

counter_1 = 0
for(s in 1:length(test))
{
  if(grepl("|", test[s]))
  {
   counter_1 = counter_1 + 1 
  }
}


#####No of complexes with 0 exp at all time points

expr_files <- list.files(path = "bnstruct/", pattern = "*expr.txt", full.names = T)
counter_expr0 <- hetero <- 0

for(i in 1:length(expr_files))
{
  expr_file <- read.table(file=expr_files[i], sep="\t", stringsAsFactors = F)
  expr_file <- sapply(expr_file[-1,-1], as.numeric)
  
  if(rowSums(expr_file)[nrow(expr_file)] == 0)
  {
    counter_expr0 = counter_expr0 + 1
    if(nrow(expr_file) > 2)
    {
      hetero = hetero + 1
    }
  }
}

