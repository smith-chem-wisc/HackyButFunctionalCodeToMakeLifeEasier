using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FusionPeptideMassVerification
{
    class PSM
    {
        public string fileName { get; set; }
        public int scan { get; set; }
        public double precursorMass { get; set; }
        public string proteinAccession { get; set; }
        public string baseSeq { get; set; }
        public string fullSeq { get; set; }
        public double[][] matchedIons { get; set; }
        public int[] matchedIonCount { get; set; }
        public double precursorMassErrorPpm { get; set; }
        public double[][] productMassErrorDa { get; set; }
        public double[][] productMassErrorPpm { get; set; }
        public double score { get; set; }
        public Boolean target { get; set; }
        public int cTarget { get; set; }
        public int cDecoy { get; set; }
        public double qValue { get; set; }
        public string nTermFrag { get; set; }
        public string cTermFrag { get; set; }
        public double expMass { get; set; }
        public double[] peakHits { get; set; }
        public int junctionIndex { get; set; }
        public List<string> bSequences { get; set; }
        public List<string> ySequences { get; set; }
        public bool NFound { get; set; }
        public bool CFound { get; set; }
        public bool isTarget { get; set; }
        public string line { get; set; }
        public int[] TD { get; set; }

        public PSM()
        {
            scan = -1;
            baseSeq = "Z";
            fullSeq = "Z";
            target = false;
            qValue = 1;
            nTermFrag = "Z";
            cTermFrag = "Z";
            junctionIndex = -1;
            expMass = -1;
            NFound = false;
            CFound = false;
        }
        public static List<PSM> parseAllPSMs(string[] candidateRead, bool keepOnlyConfidentIDs)
        {
            string[] header = candidateRead[0].Split('\t').ToArray(); //assume both files have identical headers

            int fileNameIndex = -1;
            int scanNumberIndex = -1;
            int scanPrecursorMassIndex = -1;
            int proteinAccessionIndex = -1;
            int baseSequenceIndex = -1;
            int matchedIonsIndex = -1;
            int matchedIonMassErrorDaIndex = -1;
            int matchedIonMassErrorPpmIndex = -1;
            int matchedIonCountsIndex = -1;
            int scoreIndex = -1;
            int precursorMassErrorPpmIndex = -1;
            int DCTIndex = -1;
            int cDecoy = -1;
            int cTarget = -1;
            int q = -1;
            for (int col = 0; col < header.Length; col++)
            {
                if (header[col].Equals("File Name"))
                    fileNameIndex = col;
                else if (header[col].Equals("Scan Number"))
                    scanNumberIndex = col;
                else if (header[col].Equals("Precursor Mass"))
                    scanPrecursorMassIndex = col;
                else if (header[col].Equals("Protein Accession"))
                    proteinAccessionIndex = col;
                else if (header[col].Equals("Base Sequence")) //"FullSequence" should be used for the detection of FPs containing PTMs and for missed cleave/nonspecific peptides containing PTMs
                    baseSequenceIndex = col;
                else if (header[col].Equals("Matched Ion Masses"))
                    matchedIonsIndex = col;
                else if (header[col].Equals("Matched Ion Mass Diff (Da)"))
                    matchedIonMassErrorDaIndex = col;
                else if (header[col].Equals("Matched Ion Mass Diff (Ppm)"))
                    matchedIonMassErrorPpmIndex = col;
                else if (header[col].Equals("Matched Ion Counts"))
                    matchedIonCountsIndex = col;
                else if (header[col].Equals("Score"))
                    scoreIndex = col;
                else if (header[col].Equals("Mass Diff (ppm)"))
                    precursorMassErrorPpmIndex = col;
                else if (header[col].Equals("Decoy/Contaminant/Target"))
                    DCTIndex = col;
                else if (header[col].Equals("Cumulative Decoy"))
                    cDecoy = col;
                else if (header[col].Equals("Cumulative Target"))
                    cTarget = col;
                else if (header[col].Equals("QValue"))
                    q = col;

            }
            List<PSM> allPsms = new List<PSM>();
            for (int i = 1; i < candidateRead.Length; i++)
            {
                string[] tempCand = candidateRead[i].Replace("\"","").Split('\t').ToArray();
                if (Convert.ToDouble(tempCand[q]) < 0.01)
                {
                    PSM tempPSM = new PSM();
                    tempPSM.fileName = tempCand[fileNameIndex];
                    tempPSM.scan = Convert.ToInt32(tempCand[scanNumberIndex]);
                    tempPSM.precursorMass = Convert.ToDouble(tempCand[scanPrecursorMassIndex]);
                    tempPSM.proteinAccession = tempCand[proteinAccessionIndex];
                    tempPSM.baseSeq = tempCand[baseSequenceIndex];
                    string matchedIons = tempCand[matchedIonsIndex].Replace("[", "").Replace("]", "");
                    string[] matchedIonArray = matchedIons.Split(';');
                    List<string> tempList1 = new List<string>();
                    foreach (string s in matchedIonArray)
                        if (!s.Equals(""))
                            tempList1.Add(s);
                    matchedIonArray = tempList1.ToArray();
                    tempPSM.matchedIons = new double[matchedIonArray.Length][];
                    for (int j = 0; j < matchedIonArray.Length; j++)
                    {
                        List<double> tempList = new List<double>();
                        string[] tempArray = matchedIonArray[j].Split(',');
                        foreach (string s in tempArray)
                        {
                            tempList.Add(Convert.ToDouble(s));
                        }
                        tempPSM.matchedIons[j] = tempList.ToArray();
                    }
                    string matchedIonsDa = tempCand[matchedIonMassErrorDaIndex].Replace("[", "").Replace("]", "");
                    string[] matchedIonErrorDaArray = matchedIonsDa.Split(';');
                    tempList1.Clear();
                    foreach (string s in matchedIonErrorDaArray)
                        if (!s.Equals(""))
                            tempList1.Add(s);
                    matchedIonErrorDaArray = tempList1.ToArray();
                    tempPSM.productMassErrorDa = new double[matchedIonErrorDaArray.Length][];
                    for (int j = 0; j < matchedIonErrorDaArray.Length; j++)
                    {
                        List<double> tempList = new List<double>();
                        string[] tempArray = matchedIonErrorDaArray[j].Split(',');
                        foreach (string s in tempArray)
                        {
                            tempList.Add(Convert.ToDouble(s));
                        }
                        tempPSM.productMassErrorDa[j] = tempList.ToArray();
                    }
                    string matchedIonsPpm = tempCand[matchedIonMassErrorPpmIndex].Replace("[", "").Replace("]", "");
                    string[] matchedIonErrorPpmArray = matchedIonsPpm.Split(';');
                    tempList1.Clear();
                    foreach (string s in matchedIonErrorPpmArray)
                        if (!s.Equals(""))
                            tempList1.Add(s);
                    matchedIonErrorPpmArray = tempList1.ToArray();
                    tempPSM.productMassErrorPpm = new double[matchedIonErrorPpmArray.Length][];
                    for (int j = 0; j < matchedIonErrorPpmArray.Length; j++)
                    {
                        List<double> tempList = new List<double>();
                        string[] tempArray = matchedIonErrorPpmArray[j].Split(',');
                        foreach (string s in tempArray)
                        {
                            if (!s.Contains("Infinity"))
                            {
                                tempList.Add(Convert.ToDouble(s));
                            }
                            else
                            {
                                tempList.Add(10000000000000000);
                            }
                        }
                        tempPSM.productMassErrorPpm[j] = tempList.ToArray();
                    }
                    string matchedIonsCount = tempCand[matchedIonCountsIndex].Replace("[", "").Replace("]", "");
                    string[] matchedIonCountArray = matchedIonsCount.Split(';');
                    tempList1.Clear();
                    foreach (string s in matchedIonCountArray)
                        if (!s.Equals(""))
                            tempList1.Add(s);
                    matchedIonCountArray = tempList1.ToArray();
                    tempPSM.matchedIonCount = new int[matchedIonCountArray.Length];
                    for (int j = 0; j < matchedIonCountArray.Length; j++)
                    {
                        tempPSM.matchedIonCount[j] = Convert.ToInt32(matchedIonCountArray[j]);
                    }
                    tempPSM.precursorMassErrorPpm = Convert.ToDouble(tempCand[precursorMassErrorPpmIndex].Split(' ')[0]);
                    tempPSM.score = Convert.ToDouble(tempCand[scoreIndex]);
                    tempPSM.target = tempCand[DCTIndex].Equals("T") ? true : false;
                    tempPSM.cTarget = Convert.ToInt32(tempCand[cTarget]);
                    tempPSM.cDecoy = Convert.ToInt32(tempCand[cDecoy]);
                    tempPSM.qValue = Convert.ToDouble(tempCand[q]);
                    allPsms.Add(tempPSM);
                }
            }
            return allPsms;
        }

        public PSM(string line, bool isTarget, string baseSequence, int[] TD)
        {
            this.line = line;
            this.isTarget = isTarget;
            this.baseSeq = baseSequence;
            this.TD = new int[] { TD[0], TD[1] };
        }

        public void setScan(int scan)
        {
            this.scan = scan;
        }
        public void setBaseSeq(string seq)
        {
            this.baseSeq = seq;
        }
        public void setFullSeq(string seq)
        {
            this.fullSeq = seq;
        }
        public void setTarget(Boolean target)
        {
            this.target = target;
        }
        public void setQValue(double q)
        {
            this.qValue = q;
        }
        public void setExpMass(double expMass)
        {
            this.expMass = expMass;
        }
        public void setPeakHits(double[] peakHits)
        {
            this.peakHits = peakHits;
        }
        public void generateFragments()
        {
            Random rnd = new Random(this.scan);
            this.junctionIndex = rnd.Next(1, baseSeq.Length); //A|BCDEFG|H
            //this.junctionIndex = 4;
            this.nTermFrag = this.baseSeq.Substring(0, this.junctionIndex);
            this.cTermFrag = this.baseSeq.Substring(this.junctionIndex, this.baseSeq.Length - this.junctionIndex);
        }

        public int getScan()
        {
            return this.scan;
        }
        public string getBaseSeq()
        {
            return this.baseSeq;
        }
        public string getFullSeq()
        {
            return this.fullSeq;
        }
        public Boolean getTarget()
        {
            return this.target;
        }
        public double getQValue()
        {
            return this.qValue;
        }
        public int getJunctionIndex()
        {
            return junctionIndex;
        }
        public string getNTerm()
        {
            return nTermFrag;
        }
        public string getCTerm()
        {
            return cTermFrag;
        }
        public double getExpMass()
        {
            return this.expMass;
        }
        public double[] getPeakHits()
        {
            return this.peakHits;
        }
        
        public void setPotentialSequences_B(string bSeq, int[] peaksHit)
        {
            
        }
    }
}
