using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace FusionPeptideMassVerification
{
    class FlashLFQPSM
    {
        public string fullSequence { get; set; }
        public string protein { get; set; }
        public List<double> intensities { get; set; }
        public int numKTotal { get; set; }
        public int numKLight { get; set; }
        public int numKHeavy { get; set; }
        public string baseSequence { get; set; }
        public int numOxidation { get; set; }

        public FlashLFQPSM(string fullSequence, string protein, List<double> intensities)
        {
            this.fullSequence = fullSequence;
            this.protein = protein;
            this.intensities = intensities;
            numKHeavy = fullSequence.Length - fullSequence.Replace("K[", "K").Length;
            numKTotal = fullSequence.Length - fullSequence.Replace("K", "").Length;
            numKLight = numKTotal - numKHeavy;
            baseSequence = fullSequence.Replace("[TandemMassTag:TMT_tag_lysine]", "");
            baseSequence = baseSequence.Replace("[Common Fixed:Carbamidomethyl of C]", "");
            baseSequence = baseSequence.Replace("[Common Variable:Oxidation of M]", "");
            numOxidation = fullSequence.Length - fullSequence.Replace("Ox", "O").Length;
        }
    }
}
