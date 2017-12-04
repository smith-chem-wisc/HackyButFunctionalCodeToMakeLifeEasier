using System;
using System.Globalization;
using System.IO;
using System.Data;
using System.Linq;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using System.Diagnostics;
using MathNet.Numerics;

namespace FusionPeptideMassVerification
{

    public partial class Form1 : Form
    {
        public static double[] MONOISOTOPIC_AMINO_ACID_MASSES;
        public static double[] AVERAGE_AMINO_ACID_MASSES;
        public static DataTable MassCandidates = new DataTable();
        public static DataTable DigestCandidates = new DataTable();
        public static DataTable ToleranceVsHits = new DataTable();
        public static DataTable FASTA = new DataTable();
        public static List<Protein> proteinFASTA = new List<Protein>();
        public static DataTable ModificationsDT = new DataTable();
        public static DataTable novorResults = new DataTable();
        public static double PrecursorMassTolerancePpm = 10; //(ppm)
        //public static double PrecursorMassToleranceDa = .02; //(Da)
        public static int IonsUsedMassVer = 2;
        public readonly int IonsUsedDigFilter = 6;
        //public static int troubleshooter = 0;
        //public static int RazorPeptidesRequiredForIDofPTM = 6;
        System.Windows.Forms.OpenFileDialog ofdBIons = new OpenFileDialog();
        System.Windows.Forms.OpenFileDialog ofdYIons = new OpenFileDialog();
        System.Windows.Forms.OpenFileDialog ofdFASTA = new OpenFileDialog();
        public static string BFileName;
        public static string YFileName;
        public static string FASTAFileName;
        public static double fixedModMass = 0.0; //carbamidomethyl is defaulted at 0.
        public static char[] AANames = new char[20] { 'A', 'D', 'E', 'G', 'F', 'L', 'S', 'Y', 'C', 'W', 'P', 'H', 'Q', 'R', 'I', 'M', 'T', 'N', 'K', 'V' }; //20 common AA


        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            //CalculateLocalFDR();
            //massSimilarity();
            //ConcatSeq();
            //DigestionVerification();
            //TestSoftware(); //use to generate fake fusions
            //Grapher();
            //findNumMatches();
            //SSRCalcManipulation();
            //PrintFusionCandidates();
            //DeNovoIdentifications();
            //TDHistogram();
            //Histogram();
            //SpectraZeroes(); //used for visualizing MS data in excel
            //AminoAcidMasses();
            //MassVerification();
            //MitoCartaParcer();
            //PTMsInBottomUpSamples();
            //FASTADecoyGenerator();
            //MassOfIntactProteins();
            //AminoAcidMasses();
            //ReadLiepeOutput();
            //ElucidateProductMassErrors();
            //InvestigatePSMCountDifferences();
            //FixMSGFMzidToTsvBug();
            //GenerateNonSpecificPeptidesFASTA();
            //CleanSpikes();
            //TestGamma();
            //ParseXTandem();
            //CompareSpectra();
            SILACRatios();
            //SeparateAggregatedFiles();
            //MessageBox.Show("SGALDVLQMKEEDVLK " + (MonoIsoptopicMass("SGALDVLQMKEEDVLK")+42.01056).ToString());
            //MessageBox.Show("PEPTIDE " + MonoIsoptopicMass("PEPTIDE").ToString());
            //MessageBox.Show("SEQUENCE " + MonoIsoptopicMass("SEQUENCE").ToString());
            //MessageBox.Show("PYLVSNVIELLDVDPNDQEEDGANIDLDSQR " + MonoIsoptopicMass("PYLVSNVIELLDVDPNDQEEDGANIDLDSQR").ToString());
            //MessageBox.Show("YLVSNVIELLDVDPNDQEEDGANIDLDSQR " + MonoIsoptopicMass("YLVSNVIELLDVDPNDQEEDGANIDLDSQR").ToString());
            //MessageBox.Show("LVSNVIELLDVDPNDQEEDGANIDLDSQR " + MonoIsoptopicMass("LVSNVIELLDVDPNDQEEDGANIDLDSQR").ToString());
            //MessageBox.Show("VSNVIELLDVDPNDQEEDGANIDLDSQR " + MonoIsoptopicMass("VSNVIELLDVDPNDQEEDGANIDLDSQR").ToString());
        }

        private void SILACRatios()
        {

            string[] CandidateRead = (System.IO.File.ReadAllLines(@"D:\170710_Desktop\Chemistry\Smith Research\MusSILAC\aggregateQuantifiedPeptidesByFullSeq_1mm.tsv"));
            List<FlashLFQPSM> psms = new List<FlashLFQPSM>();
            for (int i = 1; i < CandidateRead.Length; i++)
            {
                string[] array = CandidateRead[i].Split('\t').ToArray();
                List<double> intensities = new List<double>();
                for (int j = 2; j < (array.Length) / 2 + 1; j++)
                    intensities.Add(Convert.ToDouble(array[j]));
                psms.Add(new FlashLFQPSM(array[0], intensities));
            }
            string[] header = CandidateRead[0].Split('\t').ToArray();
            int numFiles = header.Length / 2 - 1;
            string[] fileNames = new string[numFiles];
            for (int k = 2; k < 2 + numFiles; k++)
                fileNames[k - 2] = header[k].Replace("Intensity_", "");
            int[] numLightOnly = new int[numFiles];
            int[] numHeavyOnly = new int[numFiles];
            List<double>[] ratios = new List<double>[numFiles];
            for (int k = 0; k < numFiles; k++)
                ratios[k] = new List<double>();
            psms = psms.OrderBy(o => o.numKHeavy).ToList();
            psms = psms.OrderBy(o => o.numOxidation).ToList();
            psms = psms.OrderBy(o => o.baseSequence).ToList();
            for (int i = 0; i < psms.Count; i++)
            {
                if(i==psms.Count-1)
                {
                    FlashLFQPSM psmLast = psms[i];
                    for (int k = 0; k < numFiles; k++)
                        if (psmLast.intensities[k] > 0)
                        {
                            if (psmLast.numKHeavy > 0)
                                numHeavyOnly[k]++;
                            if (psmLast.numKLight > 0)
                                numLightOnly[k]++;
                        }
                    break;
                }
                FlashLFQPSM psm = psms[i];
                FlashLFQPSM otherPsm = psms[i + 1];
                //if no Lysine
                if (psm.numKTotal == 0)
                    continue;

                //if not a ratio
                if (!psm.baseSequence.Equals(psms[i + 1].baseSequence) || psm.numOxidation != psms[i + 1].numOxidation)
                {
                    for (int k = 0; k < numFiles; k++)
                        if (psm.intensities[k] > 0)
                        {
                            if (psm.numKHeavy > 0)
                                numHeavyOnly[k]++;
                            else
                                numLightOnly[k]++;
                        }
                }
                else //ratio
                {
                    if (psm.numKHeavy != 0 && otherPsm.numKHeavy != 0)
                        for (int k = 0; k < numFiles; k++)
                            otherPsm.intensities[k] += psm.intensities[k];
                    else
                    {
                        for (int k = 0; k < numFiles; k++)
                        {
                            if (psm.intensities[k] > 0)
                            {
                                if (otherPsm.intensities[k] > 0)
                                    ratios[k].Add(psm.intensities[k] / (psm.intensities[k] + otherPsm.intensities[k]));
                                else
                                    numLightOnly[k]++;
                            }
                            else if (otherPsm.intensities[k] > 0)
                                numHeavyOnly[k]++;
                        }
                        if (!(i + 2 < psms.Count && otherPsm.baseSequence.Equals(psms[i + 2].baseSequence) && otherPsm.numOxidation == psms[i + 2].numOxidation))
                            i++;
                    }
                }
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\170710_Desktop\Chemistry\Smith Research\MusSILAC\SILAC_Output.txt"))
            {
                file.WriteLine("FileName" + '\t' + "Number of Peptides Observed Without a Heavy Observation" + '\t' + "Number of Peptides Observed Without a Light Observation" + '\t' + "Number of Peptides Observed as Both" + '\t' + "Average Intensity Ratio" + '\t' + "Standard Devation"+'\t'+"Percent Heavy Incorporation");
                for (int k = 0; k < numFiles; k++)
                {
                    double sum = 0;
                    foreach (double d in ratios[k])
                        sum += d;
                    double average = sum / ratios[k].Count;
                    double standardDev = 0;
                    foreach (double d in ratios[k])
                        standardDev += (average - d) * (average - d);
                    standardDev = standardDev / ratios[k].Count;
                    standardDev = Math.Pow(standardDev, 0.5);
                    file.WriteLine(fileNames[k] + '\t' + numLightOnly[k] + '\t' + numHeavyOnly[k] + '\t' + ratios[k].Count + '\t' + average + '\t' + standardDev+'\t'+(1-average)*100);
                }
            }

        }

        private void CompareSpectra()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"D:\170710_Desktop\Chemistry\Smith Research\CompressSpectra\Comparison.txt"));
            List<Tuple<double, double>> one = new List<Tuple<double, double>>();
            List<Tuple<double, double>> two = new List<Tuple<double, double>>();
            foreach (string s in CandidateRead)
            {
                string[] array = s.Split('\t').ToArray();
                if (array[0].Length > 0)
                    one.Add(new Tuple<double, double>(Convert.ToDouble(array[0]), Convert.ToDouble(array[1])));
                if (array.Length > 4 && array[3].Length > 0)
                    two.Add(new Tuple<double, double>(Convert.ToDouble(array[3]), Convert.ToDouble(array[4])));
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\170710_Desktop\Chemistry\Smith Research\CompressSpectra\commonPeaks.txt"))
            {
                int i = 0;
                int j = 0;
                while(i<one.Count && j<two.Count)
                {
                    double t = one[i].Item1 * (1 - 20.0 / 1000000.0);
                    double tt = two[j].Item1;
                    double ttt = one[i].Item1 * (1 + 20.0 / 1000000.0);
                    if (one[i].Item1 * (1 - 20.0 / 1000000.0) < two[j].Item1 && one[i].Item1 * (1 + 20.0 / 1000000.0) > two[j].Item1)
                    {
                        file.WriteLine(((one[i].Item1 + two[j].Item1) / 2 ).ToString()+ '\t' + (one[i].Item2).ToString()+'\t' + two[j].Item2);
                        i++;
                        j++;
                    }
                    else if (one[i].Item1 > two[j].Item1)
                        j++;
                    else
                        i++;
                }
            }
        }

        private void ParseXTandem()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"D:\170710_Desktop\Chemistry\Smith Research\20130504_EXQ3_MiBa_SA_Fib-2.t.xml"));
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\170710_Desktop\Chemistry\Smith Research\20130504_EXQ3_MiBa_SA_Fib-2.t.txt"))
            {
                file.WriteLine("scan" + '\t' + "experimental precursor mass" + '\t' + "charge" + '\t' + "retention time" + '\t'+ "e-value" + '\t' +"name"+'\t' + "theoretical precursor mass" + '\t' + "hyperscore" + '\t' + "sequence");
                bool nextPSM = true;
                bool catchTwo = true;
                string builder = "";
                foreach (string s in CandidateRead)
                {
                    if (s.Length < 10)
                        continue;
                    if (nextPSM)
                    {
                        if (s.Substring(0, 9).Equals("<group id"))
                        {                           
                            string[] splitString = s.Split(' ').ToArray();
                            List<int[]> indexes = new List<int[]> { new int[] { 1, 4 },new int[] { 2, 4 }, new int[] { 3, 3 },new int[] { 4, 4 },new int[] { 5, 8 } , new int[] { 6, 7 } };
                            foreach (int[] i in indexes)
                            {
                                builder += splitString[i[0]].Substring(i[1], splitString[i[0]].Length - i[1]-1);
                                builder += '\t';
                            }
                            nextPSM = false;
                        }
                    }
                    else
                    {
                        if (catchTwo && s.Substring(0, 10).Equals("<domain id"))
                        {
                            string[] splitString = s.Split(' ').ToArray();
                            List<int[]> indexes = new List<int[]> { new int[] { 5, 4 }, new int[] { 7, 12 }, new int[] { 15, 5 } };
                            foreach (int[] i in indexes)
                            {
                                builder += splitString[i[0]].Substring(i[1], splitString[i[0]].Length - i[1] - 1);
                                builder += '\t';
                            }
                            catchTwo = false;
                        }
                        else if(s.Length>25&&s.Substring(0,26).Equals("<note label=\"Description\">"))
                        {
                            builder += s.Split('.').ToArray()[1];
                            file.WriteLine(builder);
                            builder = "";
                            catchTwo = true;
                            nextPSM = true;
                        }
                    }
                }
            }
        }

        private void TestGamma()
        {
            double count = 1000000;
            double product = 5000000;
            double                    maximumLikelihood= (1.0d / count * product);

            double prevalue = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, 15);
            double prevaluea = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, 5);
            double prevalueb = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, 1);
            double prevaluec = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, 0.5);
            double prevalued = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, 0);
            double four = (1 - Math.Pow(prevalue, (count)));
            Decimal eValue = (Convert.ToDecimal((count) * four));
            double fourp = Math.Pow(count, prevalue);
            Decimal eValuep = Convert.ToDecimal(1 - (fourp / count));
            double foura = Math.Pow(count, prevaluea);
            Decimal eValuea = Convert.ToDecimal(1 - (foura / count));
            double fourb = Math.Pow(count, prevalueb);
            Decimal eValueb = Convert.ToDecimal(1 - (fourb / count));
            double fourc = Math.Pow(count, prevaluec);
            Decimal eValuec = Convert.ToDecimal(1 - (fourc / count));
            double fourd = Math.Pow(count, prevalued);
            Decimal eValued = Convert.ToDecimal(1 - (fourd / count));

        }

        private void CleanSpikes()
        {
            LoadFASTA(@"D:\170710_Desktop\Chemistry\Smith Research\Fusion Peptides\171026_MHC-I\Heck\171110_ClassicSearchOutputHeck.fasta");
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"D:\170710_Desktop\Chemistry\Smith Research\Fusion Peptides\171026_MHC-I\Heck\171110_ClassicSearchOutputSpikesHeck.txt"));
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\170710_Desktop\Chemistry\Smith Research\Fusion Peptides\171026_MHC-I\Heck\171110_ClassicSearchOutputSpikesHeckCleaned.txt"))
            {
                foreach (string s in CandidateRead)
                {
                    string sequence = s.Split('\t').ToArray()[1];
                    if (!proteinFASTA.AsParallel().Any(x => x.sequence.Contains(sequence)))
                        file.WriteLine(s);
                }
            }
        }

        private static void SeparateAggregatedFiles()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"D:\170710_Desktop\Chemistry\Smith Research\Vocal Fold\MaxQuantResults\msms07.tsv"));
            List<List<string>> lines = new List<List<string>>();
            List<string> files = new List<string>();
            for(int i=1; i<CandidateRead.Length; i++)
            {
                string[] array = CandidateRead[i].Split('\t').ToArray();
                bool done = false;
                for(int j=0; j<files.Count(); j++)
                {
                    if(files[j].Equals(array[0]))
                    {
                        lines[j].Add(CandidateRead[i]);
                        done = true;
                        break;
                    }
                }
                if (!done)
                {
                    files.Add(array[0]);
                    lines.Add(new List<string> { CandidateRead[i] });
                }
            }
            for(int j=0; j<files.Count; j++)
            {
                string s = files[j];
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\170710_Desktop\Chemistry\Smith Research\Vocal Fold\MaxQuantResults\" + s + ".tsv"))
                {
                    file.WriteLine(CandidateRead[0]);
                    foreach (string ss in lines[j])
                        file.WriteLine(ss);
                }

            }
        }

        /*   private static Dictionary<double, char> massesToResidues = new Dictionary<double, char>();

                       for(int i=0; i<Residue.ResidueMonoisotopicMass.Length; i++)
               {
                   if(!double.IsNaN(Residue.ResidueMonoisotopicMass[i]))
                   {
                       massesToResidues.Add(Residue.ResidueMonoisotopicMass[i], (char) i);
                   }
   }
   private void testMaxValue()
           {
               List<int> testints = new List<int> { 1, 2, 3, 5, 7, 6, 4, 3, 1 };
               testints.Remove(1);
               List<int> testing = new List<int>(20);
               int help = testing.Count;
               help = testing.Count();
               int asdf = 1;

               byte[] testarray = new byte[10000000];
               for (int i = 0; i < testarray.Length; i++)
                   for(int j=0; j<i%255; j++)
                       testarray[i]++;
               List<int> test = new List<int>();
               for (int i = 0; i < testarray.Length; i++)
               {
                   test.Add(i);
                   test.Clear();
               }

               int hello = 0;

               foreach(int i in test)
                   if(testarray[i]<255)
                   { }

               int highestScore = 0;

           }*/

        private void GenerateNonSpecificPeptidesFASTA()
        {
            LoadFASTA(@"\\bison\share\users\Zach\EPIC_Paper\MHC-II\170907_Human_Canonical_Uniprot.fasta");

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"\\bison\share\users\Zach\EPIC_Paper\MHC-II\170907_Human_Canonical_Uniprot_Digested_7-25_For_Rev2.fasta"))
            {
                string accession = "0";
                int total = FASTA.Rows.Count;
                int temp = 0;
                foreach (DataRow row in FASTA.Rows)
                {
                    temp++;
                    string protein = row[1].ToString();
                    //string accession = row[0].ToString();
                    int minLength = 7;
                    int maxLength = 25;
                    for (int length = minLength; length <= maxLength; length++)
                    {
                        for (int index = 0; index + length < protein.Length; index++)
                        {
                            file.WriteLine(">sp|" + accession + "|");
                            file.WriteLine(protein.Substring(index, length));

                            char[] temparray = protein.Substring(index, length).ToCharArray();
                            Array.Reverse(temparray);
                            file.WriteLine(">sp|" + accession + "R|");
                            file.WriteLine(new string(temparray));
                            accession = UpdateAccession(accession);
                        }
                    }
                }
            }
        }

        private string UpdateAccession(string accession)
        {
            char[] array = accession.ToCharArray();
            for (int i = array.Length-1; i >=0; i--)
            {
                array[i]++;
                if (array[i] == 58)
                    array[i] = 'A';
                else if (array[i] == 82) //skipR
                    array[i]++;
                else if (array[i] == 91)
                    array[i] = 'a';
                if (array[i] == 123)
                    array[i] = '0';
                else
                    return new string(array);
            }
            return "0" + new string(array);
        }

        private void FixMSGFMzidToTsvBug() //9/15/17
        {
            int i = 0;
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"\\bison\share\users\Zach\EPIC_Paper\MHC-II\TimeTrials\Task1Search\e001323-Calibrated.msgf.mzid"));
            Dictionary<string, int> dict = new Dictionary<string, int>();
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"\\bison\share\users\Zach\EPIC_Paper\MHC-II\TimeTrials\Task1Search\e001323-CalibratedFixed.msgf3.mzid"))
            {
                int j = 0;
                foreach (string s in CandidateRead)
                {
                    if (s.Contains("<Peptide id="))
                    {
                        string[] sArray = s.Split('"');
                        dict.Add(sArray[1], i);
                        file.WriteLine(sArray[0] + '"' + i + '"' + sArray[2]);
                        i++;
                    }
                    else if (s.Contains("peptide_ref=") && !s.Contains("PeptideEvidence"))
                    {
                        string[] sArray = s.Split(' ');
                        foreach (string sa in sArray)
                        {
                            if (sa.Contains("peptide_ref="))
                            {
                                string[] saArray = sa.Split('"');
                                file.WriteLine(saArray[0] + '"' + dict[saArray[1]] + '"' + saArray[2]);
                            }
                            file.Write(sa + " ");
                        }
                        file.Write('\n');
                    }
                    else
                        file.WriteLine(s);
                    j++;
                }
            }
        }

        private void InvestigatePSMCountDifferences()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"\\bison\share\users\Zach\EPIC_Paper\MHC-II\TimeTrials\Copy of Comparison.txt"));
            List<List<double>[]> listOfIDs = new List<List<double>[]>();

            foreach (string line in CandidateRead)
            {
                List<double> singleC = new List<double>();
                List<double> ionCount = new List<double>();
                List<double> classicNS = new List<double>();

                if (line.Equals(CandidateRead[0]))
                    continue;
                string newLine = line;
                newLine = newLine.Replace("[];", "");
                newLine = newLine.Replace(";[]", "");
                newLine = newLine.Replace("[", "");
                newLine = newLine.Replace("];]", "");
                newLine = newLine.Replace("]", "");
                newLine = newLine.Replace(";", ",");
                newLine = newLine.Replace("\"", "");

                if (listOfIDs.Count == 6288)
                { break; }
                try
                {
                    string[] newLineArray = newLine.Split('\t');
                    newLineArray[20].Split(',').ToList().ForEach(x => singleC.Add(Convert.ToDouble(x)));
                    newLineArray[21].Split(',').ToList().ForEach(x => ionCount.Add(Convert.ToDouble(x)));
                    newLineArray[22].Split(',').ToList().ForEach(x => classicNS.Add(Convert.ToDouble(x)));

                    if (newLineArray[20].Split(',').Count() != newLineArray[22].Split(',').Count())
                    { }
                }
                catch
                {

                }
                listOfIDs.Add(new List<double>[] { singleC, ionCount, classicNS });
                
            }

            List<double> allSingleCValues = new List<double>();
            List<double> allClassicValues = new List<double>();
            foreach(List<double>[] id in listOfIDs)
            {
                foreach (double d in id[0])
                    allSingleCValues.Add(d);// Math.Abs(d));
                foreach (double d in id[2])
                    allClassicValues.Add(d);// Math.Abs(d));
            }

            //if (allSingleCValues.Count == allClassicValues.Count) //should be true
            {
                double sumSingleC = 0;
                double sumClassic = 0;
                foreach (double d in allSingleCValues)
                    sumSingleC += d;
                foreach (double d in allClassicValues)
                    sumClassic += d;
                double avgSingleC = sumSingleC / allSingleCValues.Count;
                double avgClassic = sumClassic / allSingleCValues.Count;

                double sumForstdevClassic = 0;
                foreach (double d in allClassicValues)
                    sumForstdevClassic += (avgClassic - d) * (avgClassic - d);
                double stdevClassic = Math.Sqrt(sumForstdevClassic / (allSingleCValues.Count - 1));
                double sumForstdevSingleC = 0;
                foreach (double d in allClassicValues)
                    sumForstdevSingleC += (avgSingleC - d) * (avgSingleC - d);
                double stdevSingleC = Math.Sqrt(sumForstdevSingleC / (allSingleCValues.Count - 1));

                int[] singleCDistribution = new int[200];
                for (int i = 0; i < 200; i++)
                    singleCDistribution[i] = 0;
                foreach (double d in allSingleCValues)
                {if(d<0)
                    { }
                    int index = Convert.ToInt16(Math.Round((d - (d % 0.0005)) / 0.0005)) + 100;
                    singleCDistribution[index]++;
                }

                int[] classicDistribution = new int[200];
                for (int i = 0; i < 200; i++)
                    classicDistribution[i] = 0;
                foreach (double d in allClassicValues)
                {
                    int index = Convert.ToInt16(Math.Round((d - (d % 0.0005)) / 0.0005)) + 100;
                    classicDistribution[index]++;
                }
            }
          //  else
                throw new Exception();
        }

        private void ElucidateProductMassErrors()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"D:\170710_Desktop\Chemistry\Smith Research\ProductTolerance\MHC-II\aggregatePSMs_20ppmAroundZeroCleaned.psmtsv.txt"));

            List<PSM> allPsms = PSM.parseAllPSMs(CandidateRead, true);

            for (int multiplier = 0; multiplier <= 1200; multiplier += 200)
            {
                //multiplier -= 200;
                ProductToleranceResult[] productToleranceResults = new ProductToleranceResult[10]; //Target, Decoy with obs,comp with Da, Ppm ie TOD,TOP,TCD,TCP,TargetPrecursor, DOD,DOP,DCD,DCP, decoyprecursor
                for (int i = 0; i < productToleranceResults.Length; i++)
                    productToleranceResults[i] = new ProductToleranceResult();
                List<double> testlist = new List<double>();
                foreach (PSM psm in allPsms)
                {
                    int index;
                    if (psm.qValue < 0.01)
                    {
                        if (psm.target)
                        {
                            index = 0;
                        }
                        else
                        {
                            index = 5;
                        }
                        HashSet<int>[] foundIndexes1 = new HashSet<int>[psm.matchedIons.Length];
                        for (int i = 0; i < foundIndexes1.Length; i++)
                            foundIndexes1[i] = new HashSet<int>();
                        for (int i = 0; i < psm.matchedIons.Length; i++)
                        {
                            for (int j = 0; j < psm.matchedIons[i].Length; j++)
                            {
                                if (psm.matchedIons[i][j] > multiplier && psm.matchedIons[i][j] < multiplier+200)
                                //if (psm.precursorMass-psm.matchedIons[i][j] > 0 && psm.precursorMass-psm.matchedIons[i][j] < 200)
                                {
                                    double test = psm.precursorMass - psm.matchedIons[i][j];
                                    testlist.Add(test);
                                    foundIndexes1[i].Add(j);
                                }
                            }
                        }

                        for (int i = 0; i < psm.productMassErrorDa.Length; i++)
                        {
                            for (int j = 0; j < psm.productMassErrorDa[i].Length; j++)
                            {
                                if (foundIndexes1[i].Contains(j))
                                {
                                    productToleranceResults[index].values.Add(psm.productMassErrorDa[i][j]);
                                    if (psm.productMassErrorDa[i][j] > -1.0000000001 && psm.productMassErrorDa[i][j] < -0.9999999999)
                                        break;
                                }
                            }
                        }
                        index++;
                        for (int i = 0; i < psm.productMassErrorPpm.Length; i++)
                        {
                            for (int j = 0; j < psm.productMassErrorPpm[i].Length; j++)
                            {
                                if (foundIndexes1[i].Contains(j))
                                {
                                    productToleranceResults[index].values.Add(psm.productMassErrorPpm[i][j]);
                                    //productToleranceResults[index].values.Add(psm.productMassErrorPpm[i][j]*psm.matchedIons[i][j]/(psm.precursorMass-psm.matchedIons[i][j]));
                                    if (psm.productMassErrorPpm[i][j] > -1.0000000001 && psm.productMassErrorPpm[i][j] < -0.9999999999)
                                        break;
                                }
                            }
                        }
                        index++;
                        for (int i = 0; i < psm.productMassErrorDa.Length; i++)
                        {
                            bool b = false;
                            for (int j = 0; j < psm.productMassErrorDa[i].Length; j++)
                            {
                                if (foundIndexes1[i].Contains(j))
                                {
                                    if (b)
                                        productToleranceResults[index].values.Add(psm.productMassErrorDa[i][j]);
                                    if (psm.productMassErrorDa[i][j] > -1.0000000001 && psm.productMassErrorDa[i][j] < -0.9999999999)
                                        b = true;
                                }
                            }
                        }
                        index++;
                        for (int i = 0; i < psm.productMassErrorPpm.Length; i++)
                        {
                            bool b = false;
                            for (int j = 0; j < psm.productMassErrorPpm[i].Length; j++)
                            {
                                if (foundIndexes1[i].Contains(j))
                                {
                                    if (b)
                                        productToleranceResults[index].values.Add(psm.productMassErrorPpm[i][j]);
                                    if (psm.productMassErrorPpm[i][j] > -1.0000000001 && psm.productMassErrorPpm[i][j] < -0.9999999999)
                                        b = true;
                                }
                            }
                        }
                        index++;
                        productToleranceResults[index].values.Add(psm.precursorMassErrorPpm);
                    }
                }
                foreach (ProductToleranceResult PTR in productToleranceResults)
                {
                    PTR.numValues = PTR.values.Count;
                    double sum = 0;
                    foreach (double d in PTR.values)
                        sum += d;
                    PTR.average = sum / PTR.numValues;
                    double sumForstdev = 0;
                    foreach (double d in PTR.values)
                        sumForstdev += (PTR.average - d) * (PTR.average - d);
                    PTR.stdev = Math.Sqrt(sumForstdev / (PTR.numValues)); //don't do N-1 for sample, as throws off comparison of error.
                }
            }
          //  foreach (ProductToleranceResult PTR in productToleranceResults)
            {
             //   MessageBox.Show(PTR.numValues + " " + PTR.average + " " + PTR.stdev);
            }

/*
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"\\bison\share\groups\Neoantigen_MHC_II\LocalFDRTest\new.psmtsv"))
            {
                file.WriteLine(CandidateRead[0]);
                foreach (KeyValuePair<int, Tuple<List<PSM>, int[]>> kvp in myDict)
                {
                    foreach (PSM psm in kvp.Value.Item1)
                    {
                        string[] tempArray = psm.line.Split('\t').ToArray();
                        string tempLine = "";
                        for (int i = 0; i < tempArray.Length; i++)
                        {
                            if (i == cTarget)
                                tempLine += psm.TD[0];
                            else if (i == cDecoy)
                                tempLine += psm.TD[1];
                            else if (i == q)
                                tempLine += Convert.ToDouble(psm.TD[1]) / Convert.ToDouble(psm.TD[1] + psm.TD[0]);
                            else
                                tempLine += tempArray[i];
                            tempLine += '\t';
                        }
                        file.WriteLine(tempLine);
                    }
                }
            }*/
        }

        private void ReadLiepeOutput()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"D:\170710_Desktop\Chemistry\Smith Research\Fusion Peptides\Biological splicing papers\spliced2.txt"));
            List<string> sequences = new List<string>();
            foreach (string s in CandidateRead)
            {
                string[] tempArray = s.Split(' ').ToArray();
                foreach (string ss in tempArray)
                {
                    if (ss.Length == 9)
                    {
                        bool pass = true;
                        char[] tempchar = ss.ToCharArray();
                        foreach(char c in tempchar)
                        {
                            if ((char.IsNumber(c)))
                                pass = false;
                        }
                        if(pass)
                            sequences.Add(ss);
                    }
                }
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\170710_Desktop\Chemistry\Smith Research\Fusion Peptides\Biological splicing papers\splicedOutput2.txt"))
            {
                foreach(string s in sequences)
                { file.WriteLine(s); }
            }
        }

        private void CalculateLocalFDR()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"\\bison\share\groups\Neoantigen_MHC_II\LocalFDRTest\1323-5.psmtsv"));

            string[] header = CandidateRead[0].Split('\t').ToArray(); //assume both files have identical headers

            int fileNameIndex = -1;
            int scanNumberIndex = -1;
            int scanPrecursorMassIndex = -1;
            int proteinAccessionIndex = -1;
            int baseSequenceIndex = -1;
            int matchedIonsIndex = -1;
            int matchedIonCountsIndex = -1;
            int scoreIndex = -1;
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
                else if (header[col].Equals("Matched Ion Counts"))
                    matchedIonCountsIndex = col;
                else if (header[col].Equals("Score"))
                    scoreIndex = col;
                else if (header[col].Equals("Decoy/Contaminant/Target"))
                    DCTIndex = col;
                else if (header[col].Equals("Cumulative Decoy"))
                    cDecoy = col;
                else if (header[col].Equals("Cumulative Target"))
                    cTarget = col;
                else if (header[col].Equals("QValue"))
                    q = col;

            }
            Dictionary<int, Tuple<List<PSM>, int[]>> myDict = new Dictionary<int, Tuple<List<PSM>, int[]>>();
            for (int i = 1; i < CandidateRead.Length; i++)
            {
                string[] tempArray = CandidateRead[i].Split('\t').ToArray();
                bool DCT = tempArray[DCTIndex].Equals("T") ? true : false;
                if (myDict.TryGetValue(tempArray[baseSequenceIndex].Length, out Tuple<List<PSM>, int[]> value))
                {
                    if (DCT)
                        value.Item2[0]++;
                    else
                        value.Item2[1]++;
                    value.Item1.Add(new PSM(CandidateRead[i], DCT, tempArray[baseSequenceIndex], value.Item2));
                }
                else
                {
                    int[] tempInt = new int[] { 0, 0 };
                    if (DCT)
                        tempInt[0]++;
                    else
                        tempInt[1]++;
                    myDict.Add(tempArray[baseSequenceIndex].Length, Tuple.Create(new List<PSM> { new PSM(CandidateRead[i], DCT, tempArray[baseSequenceIndex], tempInt) }, tempInt));
                }
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"\\bison\share\groups\Neoantigen_MHC_II\LocalFDRTest\new.psmtsv"))
            {
                file.WriteLine(CandidateRead[0]);
                foreach (KeyValuePair<int, Tuple<List<PSM>, int[]>> kvp in myDict)
                {
                    foreach (PSM psm in kvp.Value.Item1)
                    {
                        string[] tempArray = psm.line.Split('\t').ToArray();
                        string tempLine = "";
                        for(int i=0; i<tempArray.Length; i++)
                        {
                            if (i == cTarget)
                                tempLine += psm.TD[0];
                            else if (i == cDecoy)
                                tempLine += psm.TD[1];
                            else if (i == q)
                                tempLine += Convert.ToDouble(psm.TD[1]) / Convert.ToDouble(psm.TD[1] + psm.TD[0]);
                            else
                                tempLine += tempArray[i];
                            tempLine += '\t';
                        }
                        file.WriteLine(tempLine);
                    }
                }
            }
        }

        private void Histogram()
        {
            double binWidth = 0.02;
            //give list of sorted values (low to high), will output list of bins with number of entries in list
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"C:\Users\zrolf\Desktop\new 3.txt"));
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\zrolf\Desktop\HistOutput.txt"))
            {
                int placeholder = 0;
                for (double bin = Convert.ToDouble(CandidateRead[0]); bin <= Convert.ToDouble(CandidateRead[CandidateRead.Count() - 1]); bin += binWidth)
                {
                    int numHits = 0;
                    for(int i=placeholder; i<CandidateRead.Count(); i++)
                    {
                        double d = Convert.ToDouble(CandidateRead[i]);
                        //MessageBox.Show(0.5544444444444.ToString() + " " + Math.Round(0.554444444, 1).ToString());
                        if(d<binWidth/2+bin&&d>bin-binWidth/2)
                        {
                            numHits++;
                        }
                        else if(d>bin+binWidth/2)
                        {
                            placeholder = i;
                            i = CandidateRead.Count();
                        }
                    }
                    file.WriteLine(bin.ToString() + '\t' + numHits.ToString());
                }
            }
        }
        private void PTMsInBottomUpSamples()
        {
            //requires DataInput.txt (headers with the accession, FullSequence, and BaseSequence of the 1% FDR peptides)
            //requires ProteinInput.txt (List of accessions that you're interested in)
            //outputs PTMOutput.txt
            LoadFASTA(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\proteasome\160927_uniprot-reviewed-human-canonical.fasta");
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Proteasome\DataInput.txt"));
            int nr = CandidateRead.GetLength(0);
            int numColumns = 3;
            string[,] data = new string[nr, numColumns]; //accession, sequenceFull, sequenceNoPTMs
            string[] lr = CandidateRead[0].Split('\t').ToArray();
            //read in headers
            data[0, 0] = lr[0];
            data[0, 1] = lr[1];
            data[0, 2] = lr[2];
            for (int r = 1; r < nr; r++)
            {
                Boolean repeat = false;
                for (int s = 0; s < r; s++)
                {
                    if (CandidateRead[s] == CandidateRead[r]) //if there is an IDENTICAL PSM, don't add this repeat
                    {
                        repeat = true;
                        data[r, 0] = "";
                        data[r, 1] = "";
                        data[r, 2] = "";
                    }
                }
                if (!repeat)
                {
                    lr = CandidateRead[r].Split('\t').ToArray();
                    data[r, 0] = lr[0];
                    data[r, 1] = lr[1];
                    data[r, 2] = lr[2];
                }
            }

            string[] Proteins = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Proteasome\ProteinInput.txt")); //accessions interested in (white space allowed)
            //string[] PTMs = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Proteasome\PTMInput.txt")); //list of observed PTMs
            List<string> PTMs = new List<string>();
            DataTable output = new DataTable();
            output.Columns.Add("Accession Numbers", typeof(string));
            foreach (string protein in Proteins)
            { output.Rows.Add(protein); }


            for (int protRow = 0; protRow < Proteins.Length; protRow++) //foreach protein of interest
            {
                string PTM = "";
                string ProtSeq = "";
                //get parent protein sequence
                foreach (DataRow row in FASTA.Rows)
                {
                    if (row[0].ToString() == Proteins[protRow])
                    {
                        ProtSeq = row[1].ToString();
                    }
                }
                //prepare indexing list
                List<List<int>> foundPTMs = new List<List<int>>(); //this is specific to each protein, while PTMs is global
                foreach (string ptm in PTMs)
                {
                    List<int> temp = new List<int>();
                    foundPTMs.Add(temp);
                }
                //List<string> usedPSMs = new List<string>(); //each PSM found for this protein
                for (int dataRow = 0; dataRow < data.Length / numColumns; dataRow++) //foreach PSM
                {
                    try
                    {
                        if (data[dataRow, 0] == Proteins[protRow]) //if PSM is from prot of interest
                        {
                            //lets find out if it has a ptm first
                            int aaIndex = 0;
                            char[] protArray = data[dataRow, 1].ToCharArray(); //parce sequence
                            for (int c = 0; c < protArray.Length; c++) //foreach aa in sequence
                            {
                                if (protArray[c] == '[')//if start of PTM
                                {
                                    c++;
                                    int bracketCounter = 1;
                                    //get the whole ptm
                                    while (bracketCounter > 0)
                                    {
                                        PTM += protArray[c];
                                        c++;
                                        if (protArray[c] == ']')
                                        {
                                            bracketCounter--;
                                        }
                                        else if (protArray[c] == '[')
                                        {
                                            bracketCounter++;
                                        }
                                    }
                                    if (!PTMs.Contains(PTM)) //if new ptm, add a column
                                    {
                                        PTMs.Add(PTM);
                                        List<int> temp = new List<int>();
                                        foundPTMs.Add(temp);
                                        output.Columns.Add(PTM, typeof(int));
                                        foreach (DataRow row in output.Rows)
                                        {
                                            /*     if (row["Accession Numbers"].ToString() == Proteins[protRow])
                                                 { row[PTM] = 1; }
                                                 else*/
                                            {
                                                row[PTM] = 0; //fill all rows of new column with zeroes.
                                            }
                                        }
                                    }
                                    for (int ptm = 0; ptm < PTMs.Count(); ptm++) //foreach ptm
                                    {
                                        if (PTMs[ptm] == PTM) //find index of list
                                        {
                                            //get index of ptm, which is index of peptide plus the location of the ptm
                                            int newPTMIndex = ProtSeq.IndexOf(data[dataRow, 2]) + aaIndex;
                                            Boolean novel = true;
                                            foreach (int oldPTMIndex in foundPTMs[ptm])
                                            {
                                                if (oldPTMIndex == newPTMIndex)
                                                {
                                                    novel = false;
                                                }
                                            }
                                            if (novel)
                                            {
                                                foundPTMs[ptm].Add(newPTMIndex);
                                            }
                                        }
                                    }
                                }
                                aaIndex++;
                                PTM = "";
                            }
                        }
                    }
                    catch { MessageBox.Show("Woops! Something was wrong with the input " + data.Length.ToString() + " " + dataRow.ToString() + " " + Proteins.Length.ToString() + " " + protRow.ToString()); }
                }
                for (int ptmIndex = 0; ptmIndex < PTMs.Count(); ptmIndex++)
                {
                    output.Rows[protRow][ptmIndex + 1] = foundPTMs[ptmIndex].Count(); //+1 is because first column is for accession numbers, not PTMs
                }
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Proteasome\PTMOutput.txt"))
            {
                string line = "Proteins";
                foreach (string ptm in PTMs)
                { line = line + "\t" + ptm; }
                file.WriteLine(line);
                for (int i = 0; i < output.Rows.Count; i++)
                {
                    line = output.Rows[i]["Accession Numbers"].ToString();
                    for (int j = 0; j < PTMs.Count; j++)
                    {
                        line = line + "\t" + output.Rows[i][PTMs[j]].ToString();
                    }
                    file.WriteLine(line);
                }
            }
        }
        private void MassOfIntactProteins()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"D:\170710_Desktop\Chemistry\Smith Research\Proteasome\ProteinIDforMass2.txt"));
            //string[,] massOutput = new string[CandidateRead.Length, 2];
            List<string> massOutput = new List<string>();
            AminoAcidMasses();
            LoadFASTA(@"D:\170710_Desktop\Chemistry\Smith Research\proteasome\160927_uniprot-reviewed-human-canonical.fasta");
            HashSet<string> found = new HashSet<string>();
            for (int i = 0; i < CandidateRead.Length; i++)
            {
                if (!found.Contains(CandidateRead[i]))
                {
                    found.Add(CandidateRead[i]);
                    for (int j = 0; j < FASTA.Rows.Count; j++)
                    {
                        if (CandidateRead[i].Equals(FASTA.Rows[j][0]))
                        {
                          //  massOutput[i, 0] = CandidateRead[i];
                            massOutput.Add( MonoIsoptopicMass(FASTA.Rows[j][1].ToString()).ToString());
                        }
                    }
                }
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\170710_Desktop\Chemistry\Smith Research\Proteasome\IntactMassOutput2.txt"))
            {
                //string line = "";
                //for (int i = 0; i < massOutput.Length; i++)
                    for (int i = 0; i < massOutput.Count; i++)
                    {
                    //line = "";
                    //for (int j = 0; j < 2; j++)
                    //{
                    //    line = line + massOutput[i, j] + "\t";
                    //}
                    //file.WriteLine(line);
                    file.WriteLine(massOutput[i]);
                }
            }
        }
        private void TDHistogram()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\TDHist.txt"));
            int highscore = Convert.ToInt16(CandidateRead[1].Split('\t').ToArray()[0]);
            int[,] histOutput = new int[highscore + 1, 3]; // score, num T, num D
            histOutput[highscore, 0] = highscore;
            histOutput[highscore, 1] = 0;
            histOutput[highscore, 2] = 0;

            int score = highscore;
            for (int i = 1; i < CandidateRead.Length; i++) //keep headers
            {
                string[] row = CandidateRead[i].Split('\t').ToArray();
                if (row[0] == score.ToString())
                {
                    if (row[1] == "D")
                    {
                        histOutput[score, 2]++;
                    }
                    else
                    {
                        histOutput[score, 1]++;
                    }
                }
                else
                {
                    score--;
                    try
                    {
                        histOutput[score, 0] = score;
                        histOutput[score, 1] = 0;
                        histOutput[score, 2] = 0;
                        i--;
                    }
                    catch { MessageBox.Show(score.ToString() + " " + i.ToString()); }
                }
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\TDHistOutput.txt"))
            {
                string line = "";
                for (int i = 0; i < histOutput.Length / 3; i++)
                {
                    line = "";
                    for (int j = 0; j < 3; j++)
                    {
                        line = line + histOutput[i, j] + "\t";
                    }
                    file.WriteLine(line);
                }
            }
        }
        private void MitoCartaParcer()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Mitochondria\MitoCartaInput.txt"));
            string[] proteinList = new string[CandidateRead.Count()];
            for (int i = 0; i < CandidateRead.Count(); i++)
            {
                string[] fusionArray = CandidateRead[i].Split('|');
                proteinList[i] = fusionArray[fusionArray.Count() - 1];
            }
            //need to choose from an xml file, not a FASTA. Have to write this code eventually... can't avoid it forever.
            string[] xmlRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Mitochondria\MitoXML.txt"));
            MessageBox.Show(xmlRead.Length.ToString());
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Mitochondria\MitoXMLOutput.txt"))
            {
                int j = 0;
                int k = 0;
                for (int i = 0; i < xmlRead.Count(); i++)
                {
                    if (xmlRead[i].Length > 6 && xmlRead[i].Substring(0, 6) == "<entry")
                    {
                        Boolean hit = false;
                        j++;
                        i++;
                        foreach (string protein in proteinList)
                        {
                            if (xmlRead[i] == "<accession>" + protein + "</accession>")
                            {
                                hit = true;
                                k++;
                                i--;
                                while (xmlRead[i] != "</entry>")
                                {
                                    file.WriteLine(xmlRead[i]);
                                    i++;
                                }
                            }
                        }
                        if (!hit) { MessageBox.Show(xmlRead[i]); }
                    }
                }
                MessageBox.Show(j.ToString() + " " + k.ToString());
            }
            //  </ entry >
            //  < entry
            /*  using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Mitochondria\MitoCartaOutput.txt"))
              {
                    file.WriteLine(output);  
              }*/
        }
        private void ConcatSeq()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\333pep.txt"));
            string[] proteinList = new string[CandidateRead.Count()];
            string concat = "";
            for (int i = 0; i < CandidateRead.Count(); i++)
            {
                concat += CandidateRead[i];
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\333FusProtSeq.txt"))
            {
                for (int i = 0; i < concat.Length; i += 60) //60 used as number of AAs per line in a FASTA file.
                {
                    string line = "";
                    if ((i + 60) < concat.Length)
                    {
                        line = concat.Substring(i, 60);
                    }
                    else
                    {
                        line = concat.Substring(i, (concat.Length - i));
                    }
                    file.WriteLine(line);
                }
            }
        }
        private void FASTADecoyGenerator() //quick and dirty. Does not account for identical sequences.
        {
            LoadFASTA(@"D:\170710_Desktop\Chemistry\Smith Research\170803_human_Uniprot_canoncial_reviewed.fasta");
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\170710_Desktop\Chemistry\Smith Research\FORWARDandREVERSED.fasta"))
            {
                foreach (DataRow row in FASTA.Rows)
                {
                    file.WriteLine(">sp|" + row[0] + "|");
                    string seq = row[1].ToString();
                    for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file.
                    {
                        string line = "";
                        if ((i + 60) < seq.Length)
                        {
                            line = seq.Substring(i, 60);
                        }
                        else
                        {
                            line = seq.Substring(i, (seq.Length - i));
                        }
                        file.WriteLine(line);
                    }
                    file.WriteLine(">sp|Decoy_" + row[0] + "|");
                    char[] charArray = seq.ToCharArray();
                    Array.Reverse(charArray);
                    seq = new string(charArray);
                    for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file.
                    {
                        string line = "";
                        if ((i + 60) < seq.Length)
                        {
                            line = seq.Substring(i, 60);
                        }
                        else
                        {
                            line = seq.Substring(i, (seq.Length - i));
                        }
                        file.WriteLine(line);
                    }
                }
            }
        }
        private void SpectraZeroes()//used for visualizing MS data in excel
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\MassImports.txt"));
            double[,] massOutput = new double[CandidateRead.Length, 2];
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\MassExports.txt"))
            {
                string line = "";
                for (int i = 0; i < CandidateRead.Length; i++)
                {
                    string[] lr = CandidateRead[i].Split('\t').ToArray();
                    double[] stuff = new double[2];
                    for (int j = 0; j < 2; j++)
                    {
                        stuff[j] = Convert.ToDouble(lr[j]);
                    }
                    line = string.Format("{0}\t{1}", (stuff[0] - 0.0001).ToString(), "0");
                    file.WriteLine(line);
                    line = string.Format("{0}\t{1}", stuff[0].ToString(), stuff[1].ToString());
                    file.WriteLine(line);
                    line = string.Format("{0}\t{1}", (stuff[0] + 0.0001).ToString(), "0");
                    file.WriteLine(line);
                }
            }
        }
        private void SSRCalcManipulation()
        {
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\SSRCalcImport.txt"));
            int numRows = CandidateRead.Length;
            double[,] SSRCalc = new double[numRows, 4];
            for (int i = 0; i < numRows; i++)
            {
                string[] input = CandidateRead[i].Split('\t');
                try
                {
                    for (int j = 0; j < 4; j++)
                    {
                        SSRCalc[i, j] = Convert.ToDouble(input[j]);
                    }
                }
                catch
                {
                    for (int j = 0; j < 4; j++)
                    {
                        if (j < 2)
                        {
                            SSRCalc[i, j] = Convert.ToDouble(input[j]);
                        }
                        else
                        {
                            SSRCalc[i, j] = 0.0;
                        }
                    }
                }
            }
            double sumY = 0;
            double sumX2 = 0;
            double sumX = 0;
            double sumXY = 0;
            for (int i = 0; i < numRows; i++)
            {
                sumY += SSRCalc[i, 1];
                sumX2 += (SSRCalc[i, 0] * SSRCalc[i, 0]);
                sumX += SSRCalc[i, 0];
                sumXY += (SSRCalc[i, 0] * SSRCalc[i, 1]);
            }
            double slope = (numRows * sumXY - (sumY * sumX)) / (numRows * sumX2 - (sumX * sumX));
            double intercept = (sumY * sumX2 - sumX * sumXY) / (numRows * sumX2 - sumX * sumX);

            double min = -10; //slide along x- axis
            double max = 220;
            double binWidth = (max - min) / 10.0;
            double binShift = 5.0;
            List<List<double>> boarderOutput = new List<List<double>>();
            for (double binStart = min; binStart < max; binStart += binShift)
            {
                int numInBin = 0;
                double sumVar = 0;
                List<double> storeBin = new List<double>();
                double minIntercept = (binStart * slope + intercept) + binStart / slope;
                double maxIntercept = ((binStart + binShift) * slope + intercept) + (binStart + binShift) / slope;
                for (int i = 0; i < numRows; i++)
                {
                    double nBord = SSRCalc[i, 0] * (-1 / slope) + maxIntercept;
                    double sBord = SSRCalc[i, 0] * (-1 / slope) + minIntercept;
                    double eBord = (SSRCalc[i, 1] - maxIntercept) / (-1 / slope);
                    double wBord = (SSRCalc[i, 1] - minIntercept) / (-1 / slope);
                    if (SSRCalc[i, 0] > wBord && SSRCalc[i, 0] < eBord && SSRCalc[i, 1] > sBord && SSRCalc[i, 1] < nBord)
                    {
                        numInBin++;
                        double xPos = SSRCalc[i, 0];
                        double yPos = SSRCalc[i, 1];
                        double tempIntercept = yPos / (-1 * slope * xPos);
                        double xLinePos = (tempIntercept - intercept) / (slope + 1 / slope);
                        double yLinePos = slope * xPos + intercept;
                        double variance = Math.Sqrt((xPos - xLinePos) * (xPos - xLinePos) + (yPos - yLinePos) * (yPos - yLinePos));
                        sumVar += variance;
                    }
                }
                //find mean variance from best fit, then create cutoff at twice that
                if (numInBin > 0)
                {
                    double meanVar = sumVar / numInBin;
                    double stdDev = Math.Sqrt(meanVar) / 2;
                    List<double> newEntry = new List<double>
                    {
                        binStart + binShift / 2
                    };
                    double yMean = (binStart + binShift / 2) * slope + intercept;
                    newEntry.Add(stdDev + yMean);
                    newEntry.Add(yMean - stdDev);
                    boarderOutput.Add(newEntry);
                }
            }
            List<Boolean> inclusionOutput = new List<Boolean>();
            for (int i = 0; i < numRows; i++) //foreach fusion
            {
                //which bin is it most central to
                double xPos = SSRCalc[i, 2];
                double yPos = SSRCalc[i, 3];
                /*double tempIntercept = yPos / (-1 * slope * xPos);
                double xLinePos = (tempIntercept - intercept) / (slope + 1 / slope);
                double yLinePos = slope * xPos + intercept;
                double variance = Math.Sqrt((xPos - xLinePos) * (xPos - xLinePos) + (yPos - yLinePos) * (yPos - yLinePos)); //a2+b2=c2*/
                double distance = 9000;
                int binNum = 0;
                for (int j = 0; j < boarderOutput.Count; j++)
                {
                    if (Math.Abs(boarderOutput[j][0] - xPos) < distance) //closest x value, not necessarily y... acceptable approximation because of small y scale
                    {
                        distance = Math.Abs(boarderOutput[j][0] - xPos);
                        binNum = j;
                    }
                }
                // double xPosBoarder = boarderOutput[binNum][0];
                double yPosBoarderHigh = boarderOutput[binNum][1];
                double yPosBoarderLow = boarderOutput[binNum][2];
                /*    double tempInterceptBoarder = yPosBoarder / (-1 * slope * xPosBoarder);
                    double xLinePosBoarder = (tempInterceptBoarder - intercept) / (slope + 1 / slope);
                    double yLinePosBoarder = slope * xPosBoarder + intercept;
                    double varianceBoarder = Math.Sqrt((xPosBoarder - xLinePosBoarder) * (xPosBoarder - xLinePosBoarder) + (yPosBoarder - yLinePosBoarder) * (yPosBoarder - yLinePosBoarder));
                    MessageBox.Show(xPos.ToString()+"_"+yPos.ToString() + "   " + xPosBoarder.ToString() + "_" + yPosBoarder.ToString());
                    MessageBox.Show(variance.ToString() + " "+varianceBoarder.ToString());
                    MessageBox.Show(xPos.ToString() + " " + xLinePos.ToString() + " " + yPos.ToString() + " " + yLinePos.ToString());
                    MessageBox.Show(xPosBoarder.ToString() + " " + xLinePosBoarder.ToString() + " " + yPosBoarder.ToString() + " " + yLinePosBoarder.ToString());*/
                Boolean withinRange = (yPos < yPosBoarderHigh && yPos > yPosBoarderLow);
                if (withinRange)
                {
                    inclusionOutput.Add(true);
                }
                else
                {
                    inclusionOutput.Add(false);
                }
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\SSRCalcBoardersExports.txt"))
            {
                foreach (List<double> list in boarderOutput)
                {
                    file.WriteLine(list[0].ToString() + '\t' + list[1].ToString() + '\t' + list[2].ToString());
                }
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\SSRCalcExports.txt"))
            {
                foreach (Boolean bol in inclusionOutput)
                {
                    file.WriteLine(bol.ToString());
                }
            }
        }
        private void NodeModule()
        {
            //given a Y ion sequence hits and spectra
            //foreach amino acid mass, 
            //if the last peak of Y ion sequence+aminoacid mass== a peak
            //save AA and mass of theoretical b ion (parent-y-18)

            //if a theo peptide does not have that b ion, do not allow it to match to that spectra.
        }

        private Tuple<List<Protein>, bool> FakeLoop(List<PSM> realPSMs, Dictionary<Protein,List<int[]>> proteinToSequences, List<Protein> proteinList)
        {
            foreach (Protein protein in proteinList)
            {
                List<int[]> foundPeptides = new List<int[]>();
                foreach (PSM psm in realPSMs)
                {
                    string pepSeq = psm.getBaseSeq();
                    int startIndex = 0;
                    string subSequence = protein.sequence;
                    while (subSequence.Contains(pepSeq))
                    {
                        if (pepSeq.Equals("AQGALANIAVDKANLEI"))
                        { }
                        startIndex = subSequence.IndexOf(pepSeq);
                        foundPeptides.Add(new int[] { startIndex, startIndex + psm.getNTerm().Length, startIndex + pepSeq.Length }); //last residue is good for length, not for positions. position is index-1;
                        subSequence = subSequence.Substring(startIndex+1, subSequence.Length - startIndex-1);
                    }
                }
                //sort values
                foundPeptides = foundPeptides.OrderBy(d => d[0]).ToList();
                if (proteinToSequences.TryGetValue(protein, out var fp))
                {
                    for(int i=0; i<fp.Count(); i=0)                
                        fp.Remove(fp[i]);
                    foreach (int[] i in foundPeptides)
                        fp.Add(i);
                }
                else
                {
                    proteinToSequences.Add(protein, foundPeptides);
                }
            }

            //generate new proteins
            List<Protein> neoProteins = new List<Protein>();
            foreach (Protein protein in proteinList)
            {
                if (protein.name.Contains("P61221"))
                { }
                if (proteinToSequences.TryGetValue(protein, out var value))
                {
                    if (value.Count == 0)
                    {
                        neoProteins.Add(protein);
                    }
                    else
                    {
                        string NSeq = "";
                        string CSeq = "";
                        int nStart = 0;
                        int cStart = 0;
                        for (int i = 0; i < value.Count; i++)
                        {
                            nStart = nStart < value[i][0] ? nStart : value[i][0];
                            NSeq += protein.sequence.Substring(nStart, value[i][1] - nStart);
                            nStart = value[i][2];

                            cStart = cStart < value[i][1] ? cStart : value[i][1];
                            if(value[i][0]>cStart)
                                CSeq += protein.sequence.Substring(cStart, value[i][0] - cStart);
                            CSeq += protein.sequence.Substring(value[i][1], value[i][2] - value[i][1]);
                            cStart = value[i][2];
                        }
                        NSeq += protein.sequence.Substring(nStart, protein.sequence.Length - nStart);
                        CSeq += protein.sequence.Substring(cStart, protein.sequence.Length - cStart);
                        neoProteins.Add(new Protein(protein.name + "_N", NSeq));
                        neoProteins.Add(new Protein(protein.name + "_C", CSeq));
                    }
                }
            }

            bool found = false;
            string JPSeq = "";
            int test = 0;
            int testtwo = 0;
            foreach(PSM psm in realPSMs)
            {
                if (neoProteins.AsParallel().Any(prot => prot.sequence.Contains(psm.nTermFrag))
                    && neoProteins.AsParallel().Any(prot => prot.sequence.Contains(psm.cTermFrag)))
                {
                    found = true;
                }
                else
                {
                    JPSeq = psm.cTermFrag + JPSeq + psm.nTermFrag;
                    test++;
                }
                if (neoProteins.AsParallel().Any(prot => prot.sequence.Contains(psm.getBaseSeq())))
                    testtwo++;
            }
            foreach(Protein protein in neoProteins)
            {
                if (protein.name.Contains("P01903"))
                { }
            }
            if(test!=0)
                neoProteins.Add(new Protein("JP_" + neoProteins.Count(), JPSeq));
            return (new Tuple<List<Protein>, bool>(neoProteins, found));
        }

        private void Grapher() //this method gets sums of numbers for each column (x4) that possess a given parameter 
        {
            //hit x4, parameter (int)
            string[] CandidateRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\FigureData.txt"));
            int nr = CandidateRead.Length;
            List<FigureData> figureData = new List<FigureData>();
            for (int i = 1; i < nr; i++)
            {
                string[] temp = CandidateRead[i].Split('\t').ToArray();
                int[] tempInt = new int[5];
                for (int j = 0; j < 5; j++)
                {
                    tempInt[j] = (Convert.ToInt16(temp[j]));
                }
                FigureData tempData = new FigureData(tempInt[0], tempInt[1], tempInt[2], tempInt[3], tempInt[4]);
                figureData.Add(tempData);
            }

            figureData.OrderBy(x => x.parameter).ToList();
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\FigureOutput.txt"))
            {
                int[] sum = { figureData[0].bClassic, figureData[0].yClassic, figureData[0].bComplementary, figureData[0].yComplementary };
                int n = 1;
                for (int i = 1; i < figureData.Count(); i++)
                {
                    if (figureData[i].parameter != figureData[i - 1].parameter)
                    {
                        file.WriteLine(figureData[i - 1].parameter.ToString() + '\t' + sum[0] * 1.0 / n + '\t' + sum[1] * 1.0 / n + '\t' + sum[2] * 1.0 / n + '\t' + sum[3] * 1.0 / n + '\t' + n);
                        sum = new int[] { figureData[i].bClassic, figureData[i].yClassic, figureData[i].bComplementary, figureData[i].yComplementary };
                        n = 1;
                    }
                    else
                    {
                        sum[0] += figureData[i].bClassic;
                        sum[1] += figureData[i].yClassic;
                        sum[2] += figureData[i].bComplementary;
                        sum[3] += figureData[i].yComplementary;
                        n++;
                    }
                }
                file.WriteLine(figureData[figureData.Count - 1].parameter.ToString() + '\t' + sum[0] * 1.0 / n + '\t' + sum[1] * 1.0 / n + '\t' + sum[2] * 1.0 / n + '\t' + sum[3] * 1.0 / n + '\t' + n);
            }
        }
        private void FindNumMatches()
        {
            //pull in real hits, need scan, length, num ions, 

            //pull in B and Y files

            //specify num cleavages

            //remove scans not in hits
            //remove extra scans

            //% of random mers accounted for in db
            AminoAcidMasses();
            LoadFASTA(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\proteasome\160927_uniprot-reviewed-human-canonical.fasta");
            List<string> db = new List<string>();
            foreach (DataRow row in FASTA.Rows)
            {
                db.Add(row[1].ToString());
            }
            int mer = 4;
            int[] indexes = new int[mer];
            for (int i = 0; i < mer; i++)
            {
                indexes[i] = 0;
            }
            int[] hits = new int[mer];
            string seq = "";
            int length = 1;

            while (indexes[0] < AANames.Count())
            {
                //get new seq
                seq = "";
                for (int n = 0; n <= length; n++)
                {
                    seq += AANames[indexes[n]];
                }

                //if new seq is present
                if (db.AsParallel().Where(x => x.Contains(seq)).Any())
                {
                    hits[length]++;
                    if (mer - 1 != length) //if not last position
                    {
                        length++; //allow m to increase
                    }
                    else //don't increment length, we're happy right now!
                    {
                        if (indexes[length] < AANames.Count() - 1) //if not last aa in possible aa
                        {
                            indexes[length-1]++;
                        }
                        else //if it is, we need to go back a bit
                        {
                            indexes[length - 1] = 0;
                            indexes[length - 2]++; //could cause crashing with weird aa
                            length--;
                        }
                    }
                }
                else //if not hit, move back one
                {
                    indexes[length - 1]++;
                }
            }
            foreach(int i in hits)
            {
                MessageBox.Show(i.ToString());
            }
        }
    
        

        private void MassSimilarity()//thought experiment on possible combinations of amino acid masses
        {
            int maxPepLength = 5;
            List<double> possibleMasses = new List<double>
            {
                0
            };
            for (int i = 0; i < maxPepLength; i++)
            {
                List<double> newMasses = new List<double>();
                foreach (double mass in possibleMasses)
                {
                    foreach (char c in AANames)
                    {
                        newMasses.Add(mass + MonoIsoptopicMass(c.ToString()));
                    }
                }
                foreach (double newMass in newMasses)
                {
                    possibleMasses.Add(newMass);
                }
            }
            possibleMasses.Sort();
            int currentMass = 0;
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\randomDistribution.txt"))
            {
                file.WriteLine(0 + '\t' + 0);
                for (double i = 0; i < maxPepLength * 300; i += 0.01)
                {
                    int numPep = 0;
                    while (possibleMasses[currentMass] > i && possibleMasses[currentMass] < i + 0.01)
                    {
                        numPep++;
                        currentMass++;
                    }
                    file.WriteLine(i + '\t' + numPep);
                    file.WriteLine(i + 0.01 + '\t' + numPep);
                    file.WriteLine(i + 0.01 + '\t' + 0);
                }
            }
        }

        //for datatable, don't use
        private void FullPSMFileReadCandidates()
        {
            MassCandidates.Columns.Add("B", typeof(string));
            MassCandidates.Columns.Add("BID", typeof(string));
            MassCandidates.Columns.Add("Y", typeof(string));
            MassCandidates.Columns.Add("YID", typeof(string));
            MassCandidates.Columns.Add("ExperimentalMass", typeof(double));
            string[] BCandidateRead = (System.IO.File.ReadAllLines(BFileName));
            string[] YCandidateRead = (System.IO.File.ReadAllLines(YFileName));
            int Bnr = BCandidateRead.GetLength(0);
            int Ynr = BCandidateRead.GetLength(0);
            DataTable DataTemp = new DataTable();
            DataTemp.Columns.Add("Filename", typeof(string));
            DataTemp.Columns.Add("Spectrum Number", typeof(int));
            DataTemp.Columns.Add("Spectrum ID", typeof(string));
            DataTemp.Columns.Add("Spectrum Title", typeof(string));
            DataTemp.Columns.Add("Retention Time (minutes)", typeof(string));
            DataTemp.Columns.Add("Precursor m/z", typeof(string));
            DataTemp.Columns.Add("Precursor Intensity", typeof(string));
            DataTemp.Columns.Add("Precursor Charge", typeof(string));
            DataTemp.Columns.Add("Precursor Mass (Da)", typeof(string));
            DataTemp.Columns.Add("Experimental Peaks", typeof(string));
            DataTemp.Columns.Add("Total Intensity", typeof(string));
            DataTemp.Columns.Add("Peptide Sequence", typeof(string));
            DataTemp.Columns.Add("Base Peptide Sequence", typeof(string));
            DataTemp.Columns.Add("Protein Description", typeof(string));
            DataTemp.Columns.Add("Start Residue Number", typeof(string));
            DataTemp.Columns.Add("Stop Residue Number", typeof(string));
            DataTemp.Columns.Add("Missed Cleavages", typeof(string));
            DataTemp.Columns.Add("Theoretical Mass(Da)", typeof(string));
            DataTemp.Columns.Add("Precursor Mass Error(Da)", typeof(string));
            DataTemp.Columns.Add("Precursor Mass Error(ppm)", typeof(string));
            DataTemp.Columns.Add("Matching Products", typeof(string));
            DataTemp.Columns.Add("Total Products", typeof(string));
            DataTemp.Columns.Add("Ratio of Matching Products", typeof(string));
            DataTemp.Columns.Add("Matching Intensity", typeof(string));
            DataTemp.Columns.Add("Fraction of Intensity Matching", typeof(string));
            DataTemp.Columns.Add("Morpheus Score", typeof(string));
            DataTemp.Columns.Add("Target?", typeof(string));
            DataTemp.Columns.Add("Decoy?", typeof(string));
            DataTemp.Columns.Add("Cumulative Target", typeof(string));
            DataTemp.Columns.Add("Cumulative Decoy", typeof(string));
            DataTemp.Columns.Add("Q-Value(%)", typeof(string));


            //for (int r = 0; r < nr; r++)
            //{
            //    string[] lr = CandidateRead[r].Split('\t').ToArray();
            //    try
            //    {
            //        lr[0] = lr[0].Split('.').ToArray()[1];
            //        lr[0] = lr[0].Replace("[carbamidomethylation of C]", "");
            //        lr[0] = lr[0].Replace("UniProt: ", "");
            //        lr[0] = RemoveNestedParentheses(lr[0], false);
            //        lr[0] = lr[0].Replace("'", "");
            //        lr[2] = lr[2].Split('.').ToArray()[1];
            //        lr[2] = lr[2].Replace("[carbamidomethylation of C]", "");
            //        lr[2] = lr[2].Replace("UniProt: ", "");
            //        lr[2] = lr[2].Replace("'", "");
            //        lr[2] = RemoveNestedParentheses(lr[2], false);
            //        //remove nested parentheses
            //        MassCandidates.Rows.Add(lr[0], lr[1].Split('|').ToArray()[1], lr[2], lr[3].Split('|').ToArray()[1], lr[4]);
            //    }
            //    catch
            //    {
            //        MessageBox.Show("ID FAIL on line " + r.ToString());
            //    }
            //}
        }

        private string RemoveNestedParentheses(string lr, bool parentheses)
        {
            //Revised 3/6/2017. MetaMorpheus now uses brackets instead of parentheses for ptms
            int nested = 0;
            char[] lrCharArray = lr.ToCharArray();
            for (int i = 0; i < lrCharArray.Count(); i++)
            {
                if (lrCharArray[i] == '[') { lrCharArray[i] = '('; } //Added 3/6/17
                else if (lrCharArray[i] == ']') { lrCharArray[i] = ')'; } //Added 3/6/17
                if (lrCharArray[i] == '(')
                {
                    if (parentheses == true)
                    {
                        lrCharArray[i] = '[';
                        nested++;
                    }
                    else
                    {
                        parentheses = true;
                    }
                }
                if (lrCharArray[i] == ')')
                {
                    if (nested > 0)
                    {
                        lrCharArray[i] = ']';
                        nested--;
                    }
                    else {
                        parentheses = false; }
                }
            }
            string cleanedSequence = new string(lrCharArray);
            return cleanedSequence;
        }

        private void LoadFASTA(string fastaFileLocation)
        {
            FASTA.Columns.Add("Name", typeof(string));
            FASTA.Columns.Add("Sequence", typeof(string));
            //string[] FASTARead=(System.IO.File.ReadAllLines("@"+FASTAFileName));
            //string[] FASTARead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\proteasome\160927_uniprot-reviewed-human-canonical.fasta"));
            //string[] FASTARead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Kojak\Release1.5.4\uniprot-rat_reviewed_130504.fasta"));
            string[] FASTARead = (System.IO.File.ReadAllLines(fastaFileLocation));
            //string[] FASTARead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\161114_Mus.fasta"));
            //string[] FASTARead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Novor Accuracy\DeNovoAccuracy\160927_uniprot-reviewed-human-canonical.fasta"));
            //string[] FASTARead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionProteinXMLSpikeGenerator\FASTASpikeOutputAndFlipped.txt"));
            int nr = FASTARead.GetLength(0);
            //MessageBox.Show(nr.ToString());
            string ProteinName = "";
            string FullSequence = "";
            for (int r = 0; r < nr; r++)
            {
                string lr = FASTARead[r];
                if (lr.StartsWith(">"))
                {
                    if (FullSequence != "") //prevents first false add.
                    {
                        //   MessageBox.Show(FullSequence);
                        FASTA.Rows.Add(ProteinName, FullSequence); //load FASTA entry into datatable
                        proteinFASTA.Add(new Protein(ProteinName, FullSequence));
                    }
                    ProteinName = lr.Substring(1); //remove >
                    FullSequence = "";
                }
                else
                {
                    FullSequence += lr;
                }
            }
            FASTA.Rows.Add(ProteinName, FullSequence); //load final FASTA entry into datatable 
            proteinFASTA.Add(new Protein(ProteinName, FullSequence));
            //MessageBox.Show("FASTA In");
            for (int i = 0; i < FASTA.Rows.Count; i++)
            {
                try
                {
                    FASTA.Rows[i][0] = FASTA.Rows[i][0].ToString().Split('|').ToArray()[1];
                    proteinFASTA[i].name = proteinFASTA[i].name.Split('|').ToArray()[1];
                }
                catch { MessageBox.Show("FASTA FAIL " + FASTA.Rows[i][0].ToString()); }
            }
        }

        private void ReadInDigestCand()
        {
            DigestCandidates.Columns.Add("BSeq", typeof(string));
            DigestCandidates.Columns.Add("BID", typeof(string));
            DigestCandidates.Columns.Add("YSeq", typeof(string));
            DigestCandidates.Columns.Add("YID", typeof(string));
            DigestCandidates.Columns.Add("ExperimentalMass", typeof(double));
            DigestCandidates.Columns.Add("SequenceMatch", typeof(string));
            string[] CandRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\DigestCandidates.txt"));
            int nr = CandRead.GetLength(0);
            for (int r = 0; r < nr; r++)
            {
                string[] lr = CandRead[r].Split('\t').ToArray();
                DigestCandidates.Rows.Add(lr[0], lr[1], lr[2], lr[3], lr[4], lr[5]);
            }
            //MessageBox.Show("Cand In");
            //dataGridView1.DataSource = Candidates;
        }

    
        private void SomeComment() {
            //start common addition/deletion/sub search, extra stuff for manual viewing
            //This is the Discovery PTM portion PTMDISCOVERY
          /*  DataView dv = Variants.DefaultView;
            dv.Sort = "MassDif asc";
            DataTable SortedVariants = dv.ToTable();
            SortedVariants.Columns.Add("PSM", typeof(double));
            SortedVariants.Columns.Add("Occurances", typeof(int));
            int checkpoint = 0;
            int NumUniqueSequences = 1;
            for (int i=1; i<SortedVariants.Rows.Count; i++)
            {
                //MessageBox.Show(SortedVariants.Rows[i][0].ToString());
                //MessageBox.Show(Convert.ToDouble(SortedVariants.Rows[checkpoint][0]).ToString());
                int counter = 0;
                //if same PTM
                if (Convert.ToDouble(SortedVariants.Rows[i][0])< Convert.ToDouble(SortedVariants.Rows[i-1][0])+0.01)//+(Convert.ToDouble(SortedVariants.Rows[checkpoint][2])*PrecursorMassTolerancePpm/1000000))
                {
                    //compare species with all others containing same PTM
                    for(int j=checkpoint;j<i;j++)
                    {
                        //if already found species (but different spectra) with same PTM
                        if (Math.Abs((Convert.ToDouble(SortedVariants.Rows[i][2]) - (Convert.ToDouble(SortedVariants.Rows[j][2])))) > 3)
                        {
                            counter++;
                        }
                        if (counter==(i-checkpoint)) //if the species didn't match with any previously found species, then record it
                        {
                            NumUniqueSequences++;
                        }
                    }
                }
                else
                {
                    double MedianMassShift = Convert.ToDouble(SortedVariants.Rows[((i - checkpoint - 1) / 2) + checkpoint][0]);
                    for (int j = checkpoint; j < i; j++)
                    {
                        SortedVariants.Rows[j][4] = MedianMassShift;
                        SortedVariants.Rows[j][5] = NumUniqueSequences;
                    }
                    //if(NumUniqueSequences> RazorPeptidesRequiredForIDofPTM+1)
                    //{
                    //    double MedianMassShift = Convert.ToDouble(SortedVariants.Rows[((i-checkpoint-1)/2)+checkpoint][0]);                        
                    //    for (int j=checkpoint; j< i;j++)
                    //    {
                    //        SortedVariants.Rows[j][3] = MedianMassShift;
                    //    }
                    //}
                    //else
                    //{
                    //    for (int j = checkpoint; j < i; j++)
                    //    {
                    //        SortedVariants.Rows[j][3] = 9999;
                    //    }
                    //}
                    checkpoint = i;
                    NumUniqueSequences = 1;
                }
                if(i==SortedVariants.Rows.Count-1)
                {
                    i++;

                    double MedianMassShift = Convert.ToDouble(SortedVariants.Rows[((i - checkpoint - 1) / 2) + checkpoint][0]);
                    for (int j = checkpoint; j < i; j++)
                    {
                        SortedVariants.Rows[j][4] = MedianMassShift;
                        SortedVariants.Rows[j][5] = NumUniqueSequences;
                    }
                }
            }
            //MessageBox.Show(SortedVariants.Rows.Count.ToString());
            dataGridView1.DataSource = SortedVariants;
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\PTMMassDifList.txt"))
            {
                foreach (DataRow PTMRow in SortedVariants.Rows)
                {
                    string line = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", PTMRow[0].ToString(), PTMRow[1].ToString(), PTMRow[2].ToString(), PTMRow[3].ToString(), PTMRow[4].ToString(), PTMRow[5].ToString());
                    file.WriteLine(line);
                }
            }
            */
        }
    

        private void DeNovoIdentifications()
        {
            ImportNovorResults();

        }

        private void ImportNovorResults()
        {
            novorResults.Columns.Add("sample", typeof(string));
            novorResults.Columns.Add("scan", typeof(int));
            novorResults.Columns.Add("score", typeof(double));
            novorResults.Columns.Add("sequence", typeof(string));
            string[] novorRead = (System.IO.File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\FusionPeptideMassVerification\NovorResults.txt"));
            int nr = novorRead.GetLength(0);
            for(int i=0; i< nr; i++)
            {
                string [] row=novorRead[i].Split('\t').ToArray();
                novorResults.Rows.Add(row[0], Convert.ToInt16(row[1]), Convert.ToDouble(row[2]), row[3]);
            }
        }


        private int LengthOfPTMs(string Sequence)
        {
            int length = 0;
            bool ModificationOn = false;
            foreach (char AA in Sequence)
            {
                if (AA == '(')
                {
                    ModificationOn = true;
                }
                if (ModificationOn == true)
                {
                    length++;
                }
                if (AA == ')')
                {
                    ModificationOn = false;
                }
            }
            return length;
        }

        private double MonoIsoptopicMass(string baseSequence)
        {
            //troubleshooter++;
            //if(Convert.ToDouble(troubleshooter)==Math.Round(Convert.ToDouble(troubleshooter)/1000)*1000)
            //{
                //MessageBox.Show(troubleshooter.ToString() + " " + baseSequence);
            //}
            double monoisotopic_mass = Constants.WATER_MONOISOTOPIC_MASS;
            Boolean ModificationOn = false;
            string ModificationName = "";
            foreach (char amino_acid in baseSequence)
            {
                //MessageBox.Show(baseSequence+" "+amino_acid.ToString()+" "+ModificationOn.ToString());
                if (amino_acid == ')') //only occurs at end of mod
                {
                    ModificationOn = false;
                    //if (ModificationName == "oxidation of M") //annotated differently than uniprot
                    //if (ModificationName == "v:Oxidation of M anywhere") //annotated differently than uniprot
                        if (ModificationName == "v:Oxidation") //Added 3/6/17
                        {
                            monoisotopic_mass += 15.99491463;
                        ModificationName = "";
                    }
                    else
                    {
                        DataRow[] PTMRow = ModificationsDT.Select("Name = '" + ModificationName + "'");
                        try
                        {
                            monoisotopic_mass += Convert.ToDouble(PTMRow[0][1]);
                        }
                        catch { MessageBox.Show("PTM " + "'" + ModificationName + "'" + " could not be found"); }
                        ModificationName = "";
                    }
                }
                if (ModificationOn == true) //only occurs if "(" already found
                {
                    ModificationName += amino_acid;
                    //MessageBox.Show(ModificationName);
                }
                if (amino_acid == '(') //start collecting PTM name
                {
                    ModificationOn = true;
                }
                if (ModificationOn == false && amino_acid!=')')
                {
                    //MessageBox.Show(baseSequence + " " + amino_acid);
                    monoisotopic_mass += GetMonoisotopicMass(amino_acid,baseSequence); //something making it here after (

                }
            }
            //MessageBox.Show(monoisotopic_mass.ToString());
            return monoisotopic_mass;
        }

        private void Test() {  
            MONOISOTOPIC_AMINO_ACID_MASSES = new double['Z' - 'A'+1];
            for (int r = 0; r < MONOISOTOPIC_AMINO_ACID_MASSES.Length; r++)
            {
                //MessageBox.Show(MONOISOTOPIC_AMINO_ACID_MASSES[r].ToString());
            }
            //MessageBox.Show(MONOISOTOPIC_AMINO_ACID_MASSES.Length.ToString());
        }

        private string baseSequence;

        public int Length { get; private set; }

        public string BaseSequence
        {
            get
            {
                return baseSequence;
            }
            private set
            {
                baseSequence = value;
                Length = value.Length;
            }
        }

        public char this[int index]
        {
            get
            {
                return baseSequence[index];
            }
        }

        private static void AminoAcidMasses()
        {
            MONOISOTOPIC_AMINO_ACID_MASSES = new double['Z' - 'A' + 1]; //makes array with 26 0's
            AVERAGE_AMINO_ACID_MASSES = new double['Z' - 'A' + 1];
            for (int i = 0; i < MONOISOTOPIC_AMINO_ACID_MASSES.Length; i++)
            {
                MONOISOTOPIC_AMINO_ACID_MASSES[i] = double.NaN;
                AVERAGE_AMINO_ACID_MASSES[i] = double.NaN;
            }

            using (StreamReader amino_acids = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "amino_acids.tsv"))) //file located in Morpheus folder
            {
                amino_acids.ReadLine();

                while (amino_acids.Peek() != -1)
                {
                    string line = amino_acids.ReadLine();
                    string[] fields = line.Split('\t');

                    char one_letter_code = char.Parse(fields[0]);
                    if (!char.IsUpper(one_letter_code))
                    {
                        throw new ArgumentOutOfRangeException("Invalid amino acid abbreviation: " + one_letter_code);
                    }
                    double monoisotopic_mass = double.Parse(fields[1], CultureInfo.InvariantCulture);
                    MONOISOTOPIC_AMINO_ACID_MASSES[one_letter_code - 'A'] = monoisotopic_mass;
                    double average_mass = double.Parse(fields[2], CultureInfo.InvariantCulture);
                    AVERAGE_AMINO_ACID_MASSES[one_letter_code - 'A'] = average_mass;
                }
            }
        }

        private void ModificationMasses()
        {
            ModificationsDT.Columns.Add("Name", typeof(string));
            ModificationsDT.Columns.Add("MonoisotopicMass", typeof(double));
            using (StreamReader uniprot_mods = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "ptmlist.txt")))
            {
                string description = null;
                string feature_type = null;
                double monoisotopic_mass_shift = double.NaN;
                double average_mass_shift = double.NaN;
                while (uniprot_mods.Peek() != -1)
                {
                    string line = uniprot_mods.ReadLine();
                    if (line.Length >= 2)
                    {
                        switch (line.Substring(0, 2))
                        {
                            case "ID":
                                description = line.Substring(5);
                                description = RemoveNestedParentheses(description, true);
                                description = description.Replace("'", "");
                                break;
                            case "FT":
                                feature_type = line.Substring(5);
                                break;
                            case "TG":
                                if (feature_type == "MOD_RES")
                                {
                                    string amino_acid = line.Substring(5);
                                }
                                break;
                            case "PP":
                                if (feature_type == "MOD_RES")
                                {
                                }
                                break;
                            case "MM":
                                monoisotopic_mass_shift = double.Parse(line.Substring(5));
                                break;
                            case "MA":
                                average_mass_shift = double.Parse(line.Substring(5));
                                break;
                            case "DR":
                                if (line.Contains("PSI-MOD"))
                                {
                                }
                                break;
                            case "//":
                                ModificationsDT.Rows.Add(description,monoisotopic_mass_shift);
                                break;
                        }
                    }
                }
            }
        }

        public static double GetMonoisotopicMass(char aminoAcid, string seq)
        {
            try {
                return MONOISOTOPIC_AMINO_ACID_MASSES[aminoAcid - 'A'];
            }
            catch {
                MessageBox.Show("Invalid amino acid '"+aminoAcid+"' caught "+seq);
                return 123; }
        }

        public static double GetAverageMass(char aminoAcid)
        {
            return AVERAGE_AMINO_ACID_MASSES[aminoAcid - 'A'];
        }

        private void DataGridView1_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {

        }

        private void Button_1(object sender, EventArgs e)
        {
            ofdBIons.Filter = "All Files (*.*)|*.*";
            ofdBIons.FilterIndex = 1;
            ofdBIons.Multiselect = false;

            if (ofdBIons.ShowDialog() == DialogResult.OK)
            {
                BFileName = ofdBIons.FileName;
                BIontxtBox.Text = ofdBIons.FileName;
            }
            //MessageBox.Show("You pressed the button");
            //if (ofdBIons.ShowDialog() == DialogResult.OK)
            //{
            //    MessageBox.Show("You pressed the button");
            //    BIontxtBox.Text = ofdBIons.FileName;
            //    ofdBIons.InitialDirectory = Path.GetDirectoryName(ofdBIons.FileName);
            //    tspbProgress.Value = tspbProgress.Minimum;

        }

        private void Button_2(object sender, EventArgs e)
        {
            ofdYIons.Filter = "All Files (*.*)|*.*";
            ofdYIons.FilterIndex = 1;
            ofdYIons.Multiselect = false;

            if (ofdYIons.ShowDialog() == DialogResult.OK)
            {
                YFileName = ofdYIons.FileName;
                YIontxtBox.Text = ofdYIons.FileName;
            }
        }

        private void BIontxtBox_TextChanged(object sender, EventArgs e)
        {

        }

        private void YIontxtBox_TextChanged(object sender, EventArgs e)
        {

        }

        private void AnalyzeButton_Click(object sender, EventArgs e)
        {
            MessageBox.Show("Close me to start!");
            //DigestionVerification();
            MessageBox.Show("Done!");
        }

        private void Label3_Click(object sender, EventArgs e)
        {

        }

        private void FASTAButton_Click(object sender, EventArgs e)
        {
            ofdFASTA.Filter = "All Files (*.*)|*.*";
            ofdFASTA.FilterIndex = 1;
            ofdFASTA.Multiselect = false;

            if (ofdFASTA.ShowDialog() == DialogResult.OK)
            {
                FASTAFileName = ofdFASTA.FileName;
                FASTAtxtBox.Text = ofdFASTA.FileName;
            }
        }
    }
}
