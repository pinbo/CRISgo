package main

import (
    "fmt"
    "path/filepath"
    "bufio"
    "os"
    s "strings"
    "strconv"
    "sort"
)

func main() {
    //command line parameters: reference file, output-name, reference-gene-name, left-flanking sequence, right flanking sequence, and sgRNA sequence
    if len(os.Args) < 7 {
        fmt.Println("This program processes all the files with extension '.fastq' in the current folder.")
        fmt.Println("Please provide 6 arguments: reference file, output name, gene name in the reference file,  left flanking sequence, right flanking sequence, and sgRNA sequence")
        os.Exit(1)
    }
    refLib := os.Args[1]
    output := os.Args[2]
    geneID:= os.Args[3]
    leftSeq:= s.ToUpper(os.Args[4])
    rightSeq:= s.ToUpper(os.Args[5])
    sgRNA:= s.ToUpper(os.Args[6])

    // fmt.Println("check reverse complement")
    // fmt.Println(ReverseComplement("ATGCggccaatt"))

    refs := get_fasta(refLib)
    refseq := s.ToUpper(refs[geneID])
    // for k := range refs {
    //     fmt.Println(">", k)
    //     fmt.Println(refs[k])
    // }
    //fmt.Println(refs["CEN5A"])
    outfile, err := os.Create(output + ".csv")
    check(err)
    defer outfile.Close()
    w := bufio.NewWriter(outfile)
    w.WriteString("fastq_file,number_of_matched_reads,number_of_reads_with_intact_gRNA,%intact_gRNA,total_indel,%indel,number_of_reads_with_leftSeq,number_of_reads_with_rightSeq,nleftSeq/nrightSeq,")
    w.WriteString("#1_indel,#1_indel_count,#2_indel,#2_indel_count,#3_indel,#3_indel_count,#4_indel,#4_indel_count,#5_indel,#5_indel_count,#6_indel,#6_indel_count,")
    w.WriteString("#7_indel,#7_indel_count,#8_indel,#8_indel_count,#9_indel,#9_indel_count,#10_indel,#10_indel_count\n")

    top10seq, err := os.Create("top10_reads_" + output + ".txt")
    check(err)
    defer top10seq.Close()
    w2 := bufio.NewWriter(top10seq)
    rawseq := getWholeSeq (refseq, leftSeq, rightSeq) // intact reference

    w2.WriteString(output + "\n")
    w2.WriteString("Reference gene ID: " + geneID + "\n")
    w2.WriteString("Left flanking sequence:  " + leftSeq + "\n")
    w2.WriteString("Right flanking sequence: " + rightSeq + "\n")
    w2.WriteString("sgRNA sequence: " + sgRNA + "\n")
    w2.WriteString("Intact reference sequence:\n" + rawseq + "\n\n")

    matches, _ := filepath.Glob("*.fastq")
    var RCmatches [] string // for 2nd round search with Reverse complement of flanking sequences
    fmt.Println("fastq file number: ", len(matches))
    for _, match := range matches {
        //fmt.Println(match)
        indels, uniqseq, nmatch, nleftSeq, nrightSeq, nsgRNA, nindel := check_fastq (match, refseq, leftSeq, rightSeq, sgRNA)
        if nmatch == 0 {
            RCmatches = append(RCmatches, match)
            continue
        }
        w.WriteString(fmt.Sprintf("%s,%d,%d,%.2f,%d,%.2f,%d,%d,%.1f,",match, nmatch, nsgRNA, float64(nsgRNA)/float64(nmatch)*100, nindel, float64(nindel)/float64(nmatch)*100, nleftSeq, nrightSeq, float64(nleftSeq)/float64(nrightSeq)))
        sorted_indels := sortMapByValue(indels)
        n := 0
        for _, kk := range sorted_indels {
            if n > 9 {
                break
            }
            n += 1
            //fmt.Println("kk", kk, "kk.Key", kk.Key, "kk.Value", kk.Value)
            w.WriteString(fmt.Sprintf("%s,%d (%.1f%%),",kk.Key, kk.Value, float64(kk.Value)/float64(nmatch)*100))
        }
        w.WriteString("\n")
        // write top 10 reads
        w2.WriteString(match + "   Total matched reads: " + strconv.Itoa(nmatch) + "\n")
        sorted_uniqseq := sortMapByValue(uniqseq)
        n = 0
        for _, kk := range sorted_uniqseq {
            if n > 9 {
                break
            }
            n += 1
            w2.WriteString(fmt.Sprintf("%s  %d\n",kk.Key, kk.Value))
        }
        w2.WriteString("\n")

    }
    // 2nd round search with Reverse complement of flanking sequences
    fmt.Println("2nd round search with Reverse complement of flanking sequences")
    leftSeqRC := ReverseComplement(leftSeq)
    rightSeqRC := ReverseComplement(rightSeq)
    sgRNARC := ReverseComplement(sgRNA)
    refseqRC := ReverseComplement(refseq)
    rawseqRC := ReverseComplement(rawseq) // intact reference
    w2.WriteString("\n\n\n----------------Check Reverse Complement of flanking sequences----------------\n\n")
    w2.WriteString("Left flanking sequence:  " + leftSeqRC + "\n")
    w2.WriteString("Right flanking sequence: " + rightSeqRC + "\n")
    w2.WriteString("sgRNA sequence: " + sgRNARC + "\n")
    w2.WriteString("Intact reference sequence:\n" + rawseqRC + "\n\n")
    for _, match := range RCmatches {
        //fmt.Println(match)
        indels, uniqseq, nmatch, nleftSeq, nrightSeq, nsgRNA, nindel := check_fastq (match, refseqRC, rightSeqRC, leftSeqRC, sgRNARC)
        if nmatch == 0 {
            continue
        }
        w.WriteString(fmt.Sprintf("%s,%d,%d,%.2f,%d,%.2f,%d,%d,%.1f,",match, nmatch, nsgRNA, float64(nsgRNA)/float64(nmatch)*100, nindel, float64(nindel)/float64(nmatch)*100, nleftSeq, nrightSeq, float64(nleftSeq)/float64(nrightSeq)))
        sorted_indels := sortMapByValue(indels)
        n := 0
        for _, kk := range sorted_indels {
            if n > 9 {
                break
            }
            n += 1
            //fmt.Println("kk", kk, "kk.Key", kk.Key, "kk.Value", kk.Value)
            w.WriteString(fmt.Sprintf("%s,%d (%.1f%%),",kk.Key, kk.Value, float64(kk.Value)/float64(nmatch)*100))
        }
        w.WriteString("\n")
        // write top 10 reads
        w2.WriteString(match + "   Total matched reads: " + strconv.Itoa(nmatch) + "\n")
        sorted_uniqseq := sortMapByValue(uniqseq)
        n = 0
        for _, kk := range sorted_uniqseq {
            if n > 9 {
                break
            }
            n += 1
            w2.WriteString(fmt.Sprintf("%s  %d\n",kk.Key, kk.Value))
        }
        w2.WriteString("\n")
    }
    // final flush
    w.Flush()
    w2.Flush()
}

// check error
func check(e error) {
    if e != nil {
        panic(e)
    }
}

// process each fastq file
func check_fastq (input string, refseq string, leftSeq string, rightSeq string, sgRNA string) (map[string]int, map[string]int, int, int, int, int, int) {
    //fmt.Println("Processing fastq file", input)
    infile, err := os.Open(input)
    check(err)
    defer infile.Close()
    scanner := bufio.NewScanner(infile)
    indels := make(map[string]int) // a dictionary to record number of different indels
    uniqseq := make(map[string]int) // a dictionary to record unique sequences from different reads 
    linenumber := 0
    nleftSeq := 0
    nrightSeq := 0
    nsgRNA := 0
    nmatch := 0 // number of reads with both leftSeq and rightSeq
    nindel := 0 // number of reads with indels
    rawseq := getWholeSeq (refseq, leftSeq, rightSeq) // intact reference
    //fmt.Println("raw seq is", rawseq)
    rawlen := len(rawseq) // length of unedited sequence
    //fmt.Println("rawlen is", rawlen)
    for scanner.Scan() {
        //fmt.Println(scanner.Text())
        linenumber += 1
        if linenumber % 4 != 2 {
            continue
        } // only check the sequence line
        line := scanner.Text() // "\n" is already trimmed
        if s.Contains(line, leftSeq) && s.Contains(line, rightSeq) {
            nmatch += 1
            nleftSeq += 1
            nrightSeq += 1
            if s.Contains(line, sgRNA) {
                nsgRNA += 1
            }
            editedseq := getWholeSeq (line, leftSeq, rightSeq)
            uniqseq[editedseq] += 1
            editedlen := len(editedseq)
            if editedlen - rawlen != 0 {
                nindel += 1
            }
            indels[strconv.Itoa(editedlen - rawlen)] += 1

        } else if s.Contains(line, leftSeq) {
            nleftSeq += 1
        } else if s.Contains(line, rightSeq) {
            nrightSeq += 1
        } else {
            continue
        }
    }
    fmt.Println("Processing", input, "matches", nmatch)
    return indels, uniqseq, nmatch, nleftSeq, nrightSeq, nsgRNA, nindel
}

// get the whole sequence when giving left and right flanking sequences
func getWholeSeq (refseq string, leftSeq string, rightSeq string) string {
    if s.Contains(refseq, leftSeq) && s.Contains(refseq, rightSeq) {
        start := s.Index(refseq, leftSeq)
        end := s.Index(refseq, rightSeq) + len(rightSeq)
        return refseq[start:end]
    } else if s.Contains(refseq, leftSeq) {
        fmt.Println("Error, rightSeq is not in the refseq")
        return ""
    } else if s.Contains(refseq, rightSeq) {
        fmt.Println("Error, leftSeq is not in the refseq")
        return ""
    } else {
        fmt.Println("Error, nither flanking sequence is in the refseq")
        return ""
    }
}

// sort map keys by value
// A data structure to hold a key/value pair.
type Pair struct {
    Key string
    Value int
}
type PairList []Pair
func (p PairList) Len() int           { return len(p) }
func (p PairList) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }
func (p PairList) Less(i, j int) bool { return p[i].Value < p[j].Value }
// A function to turn a map into a PairList, then sort and return it. 
func sortMapByValue(m map[string]int) PairList {
    p := make(PairList, len(m))
    i := 0
    for k, v := range m {
        p[i] = Pair{k, v}
        i++
    }
    sort.Sort(sort.Reverse(p))
    return p
}

// read fasta file
func get_fasta (input string) map[string]string {
    infile, err := os.Open(input)
    check(err)
    defer infile.Close()
    scanner := bufio.NewScanner(infile)
    fasta := map[string]string{}
    seqName := "abc"
    for scanner.Scan() {
        //fmt.Println(scanner.Text())
        line := scanner.Text() // "\n" is already trimmed
        if s.HasPrefix(line, ">") {
            line := s.TrimLeft(line, "> ") // in case "> seq1"
            ll := s.Split(line, " ")
            seqName = ll[0]
            fasta[seqName] = ""
        } else {
            fasta[seqName] += line
        }
    }
    return fasta
}

// Reverse returns its argument string reversed rune-wise left to right.
func Reverse(s string) string {
    r := []rune(s)
    for i, j := 0, len(r)-1; i < len(r)/2; i, j = i+1, j-1 {
        r[i], r[j] = r[j], r[i]
    }
    return string(r)
}

// Reverse complement a sequence
func ReverseComplement(seq string) string {
    s1 := "BDHKMNRSVWYATGCbdhkmnrsvwyatgc"
    s2 := "VHDMKNYSBWRTACGvhdmknysbwrtacg"
    seqDict := map[string]string{}
    for i, c := range s1 {
        seqDict[string(c)] = string(s2[i])
    }
    ss := s.Split(seq, "")
    // reverse the array s
    for i, j := 0, len(ss)-1; i < j; i, j = i+1, j-1 {
        ss[i], ss[j] = ss[j], ss[i]
    }
    for i, c := range ss {
        ss[i] = seqDict[c]
    }

    return s.Join(ss, "")
}
