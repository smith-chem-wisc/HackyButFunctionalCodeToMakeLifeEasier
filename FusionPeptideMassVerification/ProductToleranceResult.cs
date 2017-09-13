using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FusionPeptideMassVerification
{
    class ProductToleranceResult
    {
        public ProductToleranceResult()
        {
            values = new List<double>();
            average = 0;
            stdev = 0;
            numValues = 0;
        }

        public List<double> values { get; set; }
        public double average { get; set; }
        public double stdev { get; set; }
        public int numValues { get; set; } 
    }
}
