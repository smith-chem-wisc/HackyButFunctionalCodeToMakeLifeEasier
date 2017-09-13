using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FusionPeptideMassVerification
{
    class FigureData
    {
        public int bClassic { get; set; }
        public int yClassic { get; set; }
        public int bComplementary { get; set; }
        public int yComplementary { get; set; }
        public int parameter { get; set; }

        public FigureData(int b, int y, int bc, int yc, int p)
        {
            bClassic = b;
            yClassic = y;
            bComplementary = bc;
            yComplementary = yc;
            parameter = p;
        }
    }
}
