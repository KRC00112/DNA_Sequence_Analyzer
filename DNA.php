<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DNA Sequence Analyzer</title>
</head>
<body>
   

<hr>
    <h1 style="font-family:'Trebuchet MS', sans-serif; font-size:30px">DNA Sequence Analyzer</h1>
<hr>
    <form action="DNA.php" method="post">
        <label style="font-family: Tahoma, sans-serif">Enter DNA Sequence: </label>    
        <input type="text" name="dna" value="<?php if(isset($_POST['dna'])) echo $_POST['dna']; ?>">
        <input type="submit" value="Submit"><br>
    </form>

<hr>

</body>
</html>


<?php

if(isset($_POST["dna"])) {
    

        //////////////////////////////////////////////////////////////Functions get declared here////////////////////////////////////////////////////////...
        ////////////////////////////////////////////////////////////function to find the nucleotide frequecy of the DNA sequence ////////////////////////////////////////////////////


    function nucleotideFrequency($seq) {
        $countA = 0;
        $countT = 0;
        $countG = 0;
        $countC = 0;
        $countU = 0;
        $nd = new SplDoublyLinkedList();
        $nucleotides = ['A', 'T', 'G', 'C', 'U'];
        for ($i = 0; $i < strlen($seq); $i++) {
            if ($seq[$i] == 'A') {
                $countA++;
            }
            if ($seq[$i] == 'T') {
                $countT++;
            }
            if ($seq[$i] == 'G') {
                $countG++;
            }
            if ($seq[$i] == 'C') {
                $countC++;
            }
            if ($seq[$i] == 'U') {
                $countU++;
            }
        }
        $nd->push("count of A: " . $countA);
        $nd->push("count of G: " . $countG);
        $nd->push("count of C: " . $countC);
        if (strpos($seq, 'U') !== false) {
            $nd->push("count of U: " . $countU);
        } else {
            $nd->push("count of T: " . $countT);
        }
        return $nd;
    }

        ////////////////////////////////////////////////////////////function to find all the codons in the DNA sequence////////////////////////////////////////////////////


    function codons($seq) {
        $seqLength = strlen($seq);
        $remainder = 0;
        if ($seqLength % 3 != 0) {
            $remainder = $seqLength % 3;
            $seqLength = $seqLength - $remainder;
        }
        $j = 3; 
        $singleCodon = "";
        $codons = array();
        for ($i = 0; $i < $seqLength; $i += 3) {
            $singleCodon = substr($seq, $i, $j); 
            $codons[] = $singleCodon;
        }
        return $codons;
    }

        ////////////////////////////////////////////////////////////function to find the frequency of codons ////////////////////////////////////////////////////


    
    function codonFrequency($codonList) {
        $stringCodons = implode(" ", $codonList);
        $syms = array("[", "]", ",");
        foreach ($syms as $i) {
            $stringCodons = str_replace($i, "", $stringCodons);
        }
        $dict = explode(" ", $stringCodons);
        $hs = array();
    
        if (empty($stringCodons)) {
            $hs[] = "empty";
            return $hs;
        } else {
            $uniqueCodons = array_unique($dict);
            foreach ($uniqueCodons as $word) {
                $count = 0;
                foreach ($dict as $codon) {
                    if (strcasecmp($codon, $word) == 0) {
                        $count++;
                    }
                }
                $hs[] = $word . ": " . $count;
            }
            return $hs;
        }
    }

        ////////////////////////////////////////////////////////////functions to find the reverse compliment of the DNA sequence////////////////////////////////////////////////////



    function reverseComplementOfDna($dnaSeq) {
        $store = "";
        $nucleotide='';
        for ($i = 0; $i < strlen($dnaSeq); $i++) {
            $nucleotide = $dnaSeq[$i];
            if ($nucleotide == 'A') {
                $nucleotide = 'T';
            } else if ($nucleotide == 'T') {
                $nucleotide = 'A';
            } else if ($nucleotide == 'G') {
                $nucleotide = 'C';
            } else if ($nucleotide == 'C') {
                $nucleotide = 'G';
            }
            $store = $store . $nucleotide;
        }
        return $store;
    }

        ////////////////////////////////////////////////////////////GC content calculation function////////////////////////////////////////////////////

    function gcContentCalculation($dnaSeq) {
        $countG = 0;
        $countC = 0;
        $countA = 0;
        $countT = 0;
        $gcContent = 0;
        for ($i = 0; $i < strlen($dnaSeq); $i++) {
            if ($dnaSeq[$i] == 'G') {
                $countG++;
            } else if ($dnaSeq[$i] == 'C') {
                $countC++;
            } else if ($dnaSeq[$i] == 'A') {
                $countA++;
            } else if ($dnaSeq[$i] == 'T') {
                $countT++;
            }
        }
        $gcContent = (($countG + $countC) / ($countG + $countC + $countT + $countA)) * 100;
        return $gcContent;
    }

    ////////////////////////////////////////////////////////////function for transcription(converts DNA to RNA)////////////////////////////////////////////////////

    function transcription($dnaSeq) {
        $rnaSeq = str_replace('T', 'U', $dnaSeq);
        return $rnaSeq;
    }


    /////////////////////////////////////////////////////////////function for translation(converts RNA to protein sequence)/////////////////////////////////////////////////////

    function translation($rnaSeq) {
        $codon = "";
        $j = 3;
        $aminoAcid = "";
        $remainder = 0;
        $rnaLength = strlen($rnaSeq);
        if ($rnaLength % 3 != 0) {
            $remainder = $rnaLength % 3;
            $rnaLength = $rnaLength - $remainder;
        }
        for ($i = 0; $i < $rnaLength; $i += 3) {
            $codon = substr($rnaSeq, $i, 3);
            if (strcasecmp($codon, "UUU") == 0 || strcasecmp($codon, "UUC") == 0) {
                $codon = "F";
            }
            if (strcasecmp($codon, "UUA") == 0 || strcasecmp($codon, "UUG") == 0 || strcasecmp($codon, "CUU") == 0 || strcasecmp($codon, "CUC") == 0 || strcasecmp($codon, "CUA") == 0 || strcasecmp($codon, "CUG") == 0) {
                $codon = "L";
            }
            if (strcasecmp($codon, "AUU") == 0 || strcasecmp($codon, "AUC") == 0 || strcasecmp($codon, "AUA") == 0) {
                $codon = "I";
            }
            if (strcasecmp($codon, "AUG") == 0) {
                $codon = "M";
            }
            if (strcasecmp($codon, "GUU") == 0 || strcasecmp($codon, "GUC") == 0 || strcasecmp($codon, "GUA") == 0 || strcasecmp($codon, "GUG") == 0) {
                $codon = "V";
            }
            if (strcasecmp($codon, "UCU") == 0 || strcasecmp($codon, "UCC") == 0 || strcasecmp($codon, "UCA") == 0 || strcasecmp($codon, "UCG") == 0 || strcasecmp($codon, "AGU") == 0 || strcasecmp($codon, "AGC") == 0) {
                $codon = "S";
            }
            if (strcasecmp($codon, "CCU") == 0 || strcasecmp($codon, "CCC") == 0 || strcasecmp($codon, "CCA") == 0 || strcasecmp($codon, "CCG") == 0) {
                $codon = "P";
            }
            if (strcasecmp($codon, "ACU") == 0 || strcasecmp($codon, "ACC") == 0 || strcasecmp($codon, "ACA") == 0 || strcasecmp($codon, "ACG") == 0) {
                $codon = "T";
            }
            if (strcasecmp($codon, "GCU") == 0 || strcasecmp($codon, "GCC") == 0 || strcasecmp($codon, "GCA") == 0 || strcasecmp($codon, "GCG") == 0) {
                $codon = "A";
            }
            if (strcasecmp($codon, "UAU") == 0 || strcasecmp($codon, "UAC") == 0) {
                $codon = "Y";
            }
            if (strcasecmp($codon, "UAA") == 0 || strcasecmp($codon, "UAG") == 0 || strcasecmp($codon, "UGA") == 0) {
                $codon = "*";
            }
            if (strcasecmp($codon, "CAU") == 0 || strcasecmp($codon, "CAC") == 0) {
                $codon = "H";
            }
            if (strcasecmp($codon, "CAA") == 0 || strcasecmp($codon, "CAG") == 0) {
                $codon = "Q";
            }
            if (strcasecmp($codon, "AAU") == 0 || strcasecmp($codon, "AAC") == 0) {
                $codon = "N";
            }
            if (strcasecmp($codon, "AAA") == 0 || strcasecmp($codon, "AAG") == 0) {
                $codon = "K";
            }
            if (strcasecmp($codon, "GAU") == 0 || strcasecmp($codon, "GAC") == 0) {
                $codon = "D";
            }
            if (strcasecmp($codon, "GAA") == 0 || strcasecmp($codon, "GAG") == 0) {
                $codon = "E";
            }
            if (strcasecmp($codon, "UGU") == 0 || strcasecmp($codon, "UGC") == 0) {
                $codon = "C";
            }
            if (strcasecmp($codon, "UGG") == 0) {
                $codon = "W";
            }
            if (strcasecmp($codon, "CGU") == 0 || strcasecmp($codon, "CGC") == 0 || strcasecmp($codon, "CGA") == 0 || strcasecmp($codon, "CGG") == 0 || strcasecmp($codon, "AGA") == 0 || strcasecmp($codon, "AGG") == 0) {
                $codon = "R";
            }
            if (strcasecmp($codon, "GGU") == 0 || strcasecmp($codon, "GGC") == 0 || strcasecmp($codon, "GGA") == 0 || strcasecmp($codon, "GGG") == 0) {
                $codon = "G";
            }
            $aminoAcid .= $codon;
            $j += 3;
        }
        return $aminoAcid;
    }

        ///////////////////////////////////////////////////////function for aminoAcid Frequency////////////////////////////////////////////////


    function aminoAcidFrequency($seq) {
        $dict = str_split($seq);
        $hs = [];
        if(empty($seq)) {
            $hs[] = "empty";
            return $hs;
        } else {
            $count = 0;
            $word = "";
            foreach($dict as $j => $word) {
                foreach ($dict as $i => $value) {
                    if (strcasecmp($value, $word) == 0) {
                        $count++;
                    }
                }
    
                $hs[$word] = $word.": ".$count;
                $count = 0;
            }
            return array_values($hs);
        }
    }

        ///////////////////////////////////////////////////////function for secondary structure affinity////////////////////////////////////////

    function secondaryStructureAffinity($proteinSeq) {
        $alphaCount = 0;
        $betaCount = 0;
        $coilCount = 0;
    
        for ($i = 0; $i < strlen($proteinSeq); $i++) {
    
            if (in_array($proteinSeq[$i], ['A', 'R', 'D', 'Q', 'E', 'G', 'I', 'L', 'K', 'M', 'F', 'P', 'W', 'V'])) {
    
                $alphaCount++;
            }
            if (in_array($proteinSeq[$i], ['N', 'C', 'Y'])) {
    
                $betaCount++;
            }
            if (in_array($proteinSeq[$i], ['H', 'S', 'T'])) {
    
                $coilCount++;
            }
            if ($proteinSeq[$i] == '*') {
    
                continue;
            }
        }
    
        $countArray = ["\u{03B1}: $alphaCount", "\u{03B2}: $betaCount", "C: $coilCount"];
        $stringArray = implode(", ", $countArray);
    
        $maxCount = max($alphaCount, $betaCount, $coilCount);
    
        $symbol = '';
    
        if ($maxCount == $betaCount) {
            $symbol = "\u{03B2}";
        } elseif ($maxCount == $alphaCount) {
            $symbol = "\u{03B1}";
        } elseif ($maxCount == $coilCount) {
            $symbol = "C";
        }
    
        return "[$stringArray]" . "[Max: $maxCount][Has affinity towards ($symbol) secondary structure]";
    }
    
    
    
    
        ///////////////////////////////////All Functions above this////////////////////////////////////////////////////////////////
        ///////////////////////////////////Call the functions below this//////////////////////////////////////////////////////////


    $dna = $_POST["dna"];

    echo"<pre>";
    if (!preg_match("/^[atgcATGC\s]+$/", $dna)) {
        echo "Enter a valid DNA Sequence";
    }
    
    
    

    else{
    
    $uppercaseDna = strtoupper($dna);
    $uppercaseDna = str_replace(" ", "", $uppercaseDna);
    $lengthOfDnaSeq=strlen($uppercaseDna);
    $nucleotideFreq=nucleotideFrequency($uppercaseDna);
    $codonsFormed=codons($uppercaseDna);
    $totalNoOfCodons=count($codonsFormed);
    $codonFreq=codonFrequency($codonsFormed);
    $reverseComp=reverseComplementOfDna($uppercaseDna);
    $gcContent = number_format(gcContentCalculation($uppercaseDna), 2);
    $rna=transcription($uppercaseDna);
    $proteinSeq=translation($rna);
    $aminoAcidFreq=aminoAcidFrequency($proteinSeq);
    $proteinSecondaryStructure=secondaryStructureAffinity($proteinSeq);

    echo "The DNA Sequence                                 : {$uppercaseDna}<br>";
    echo "The Length of the DNA Sequence                   : {$lengthOfDnaSeq} Nucleotides<br>";
    echo "The Nucleotide Frequency of the DNA Sequence     : [".implode(", ", iterator_to_array($nucleotideFreq))."]<br>";
    echo "Codons That can be formed from the DNA Sequence  : [".implode(", ", $codonsFormed)."][{$totalNoOfCodons}]<br>";
    echo "The Codon Frequency of the DNA Sequence          : [".implode(", ", $codonFreq)."]<br>";
    echo "The Reverse Complement of the DNA Sequence       : {$reverseComp}<br>";
    echo "The GC-Content of the DNA Sequence               : {$gcContent} %<br>";
    echo "The RNA Sequence                                 : {$rna}<br>";
    echo "The Protein Sequence                             : {$proteinSeq}<br>";
    echo "The Amino Acid Frequency of the Protein Sequence : [".implode(", ", $aminoAcidFreq)."]<br>";
    echo "The Analysis of the Protein Secondary Structure  : {$proteinSecondaryStructure} <br>";
    
    }
    echo "</pre>";
    
}
?>