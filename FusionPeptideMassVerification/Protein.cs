using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FusionPeptideMassVerification
{
    public class Protein
    {
        public string name { get; set; }
        public string sequence { get; set; }

        public Protein(string name, string sequence)
        {
            this.name = name;
            this.sequence = sequence;
        }
    }
}
