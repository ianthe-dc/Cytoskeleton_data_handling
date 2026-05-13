library(flextable)


antidf <- data.frame("Antibody/stain" = c("Anti-vimentin", "Anti α-tubulin", "Anti phospho-H2AX", "Alexa Fluor phalloidin 488", "Alexa Fluor phalloidin 633", "Alexa Fluor 594 anti-Rabbit IgG", "Alexa Fluor 488 anti-Mouse IgG", "Alexa Fluor 488 anti-Rat IgG"),
                     "Primary (P) or Secondary (S)" = c("P", "P", "P", "N/A", "N/A", "S", "S", "S"),
                     Host = c("Mouse", "Rat", "Rabbit", "N/A", "N/A", "Chicken", "Goat", "Donkey"),
                     Clonality = c("Monoclonal", "Monoclonal", "Monoclonal", "N/A", "N/A", "Polyclonal", "Polyclonal", "Polyclonal"),
                     Manufacturer = c("Invitrogen", "Bio-Rad", "Cell Signaling Technology", "ThermoFisher", "ThermoFisher", "Invitrogen", "Invitrogen", "Invitrogen"),
                     "Cat. number" = c("MA5-11883", "MCA78G", "9718S", "A12379", "A22284", "A21442", "A11001", "A21208"),
                     "IF conc." = c("1:100", "1:200", "1:250", "1:80", "1:80", "1:500", "1:500", "1:500"),
                     check.names = FALSE)

antidf

antibodies <- flextable(antidf)


antibodies

save_as_image(antibodies, path = "/Users/vbss76/Desktop/CytoAge/Thesis/antibodies.png", dpi = 300, width = 130, height = 140, units = "mm")



qpcr <- data.frame(Gene = c("p16 pair 11", "p21 pair 6", "RPLP0 pair 2", "VAMP7 pair 2", "L3MBTL2 pair 2"),
                   "Forward sequence" = c("5’-CCAACGCACCGAATAGTTACG-3’", "5’- TCAGTTCCTTGTGGAGCCGGA-3’", "5’-CGTCCTCGTGGAAGTGACAT-3’", "5’-TTTAATTGCTCAAAGCTGTCGCC-3’", "5’-ACACCAACAGGACAAGACGCTC-3’"),
                   "Reverse sequence" = c("5’-CGGGGATGTCTGAGGGACCTT-3’", "5’- GTTCTGACATGGCGCCTCCT-3’", "5’-CCGGATATGAGGCAGCAGTT-3’", "5’-AATGGGCTAGCAAAATGTTCCC-3’", "5’-AGCCTTCATACCGAAGCAGCAC-3’"),
                   check.names = FALSE)

primers <- flextable(qpcr)

primers

save_as_image(primers, path = "/Users/vbss76/Desktop/CytoAge/Thesis/Figures/primers.png", dpi = 300, width = 130, height = 140, units = "mm")
