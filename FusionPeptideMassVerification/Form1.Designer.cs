namespace FusionPeptideMassVerification
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.dataGridView1 = new System.Windows.Forms.DataGridView();
            this.splitContainer1 = new System.Windows.Forms.SplitContainer();
            this.label3 = new System.Windows.Forms.Label();
            this.FASTAButton = new System.Windows.Forms.Button();
            this.FASTAtxtBox = new System.Windows.Forms.TextBox();
            this.AnalyzeButton = new System.Windows.Forms.Button();
            this.YIontxtBox = new System.Windows.Forms.TextBox();
            this.BIontxtBox = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.label1 = new System.Windows.Forms.Label();
            this.YIonButton = new System.Windows.Forms.Button();
            this.BIonButton = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer1)).BeginInit();
            this.splitContainer1.Panel1.SuspendLayout();
            this.splitContainer1.Panel2.SuspendLayout();
            this.splitContainer1.SuspendLayout();
            this.SuspendLayout();
            // 
            // dataGridView1
            // 
            this.dataGridView1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.dataGridView1.Location = new System.Drawing.Point(0, 0);
            this.dataGridView1.Name = "dataGridView1";
            this.dataGridView1.RowTemplate.Height = 24;
            this.dataGridView1.Size = new System.Drawing.Size(885, 490);
            this.dataGridView1.TabIndex = 0;
            this.dataGridView1.CellContentClick += new System.Windows.Forms.DataGridViewCellEventHandler(this.DataGridView1_CellContentClick);
            // 
            // splitContainer1
            // 
            this.splitContainer1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.splitContainer1.Location = new System.Drawing.Point(0, 0);
            this.splitContainer1.Name = "splitContainer1";
            // 
            // splitContainer1.Panel1
            // 
            this.splitContainer1.Panel1.Controls.Add(this.label3);
            this.splitContainer1.Panel1.Controls.Add(this.FASTAButton);
            this.splitContainer1.Panel1.Controls.Add(this.FASTAtxtBox);
            this.splitContainer1.Panel1.Controls.Add(this.AnalyzeButton);
            this.splitContainer1.Panel1.Controls.Add(this.YIontxtBox);
            this.splitContainer1.Panel1.Controls.Add(this.BIontxtBox);
            this.splitContainer1.Panel1.Controls.Add(this.label2);
            this.splitContainer1.Panel1.Controls.Add(this.label1);
            this.splitContainer1.Panel1.Controls.Add(this.YIonButton);
            this.splitContainer1.Panel1.Controls.Add(this.BIonButton);
            // 
            // splitContainer1.Panel2
            // 
            this.splitContainer1.Panel2.Controls.Add(this.dataGridView1);
            this.splitContainer1.Size = new System.Drawing.Size(1333, 490);
            this.splitContainer1.SplitterDistance = 444;
            this.splitContainer1.TabIndex = 1;
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(21, 177);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(78, 17);
            this.label3.TabIndex = 9;
            this.label3.Text = "FASTA File";
            this.label3.Click += new System.EventHandler(this.Label3_Click);
            // 
            // FASTAButton
            // 
            this.FASTAButton.Location = new System.Drawing.Point(362, 195);
            this.FASTAButton.Name = "FASTAButton";
            this.FASTAButton.Size = new System.Drawing.Size(69, 24);
            this.FASTAButton.TabIndex = 8;
            this.FASTAButton.Text = "Browse";
            this.FASTAButton.UseVisualStyleBackColor = true;
            this.FASTAButton.Click += new System.EventHandler(this.FASTAButton_Click);
            // 
            // FASTAtxtBox
            // 
            this.FASTAtxtBox.Location = new System.Drawing.Point(24, 197);
            this.FASTAtxtBox.Name = "FASTAtxtBox";
            this.FASTAtxtBox.Size = new System.Drawing.Size(316, 22);
            this.FASTAtxtBox.TabIndex = 7;
            // 
            // AnalyzeButton
            // 
            this.AnalyzeButton.Location = new System.Drawing.Point(176, 244);
            this.AnalyzeButton.Name = "AnalyzeButton";
            this.AnalyzeButton.Size = new System.Drawing.Size(83, 28);
            this.AnalyzeButton.TabIndex = 6;
            this.AnalyzeButton.Text = "Analyze";
            this.AnalyzeButton.UseVisualStyleBackColor = true;
            this.AnalyzeButton.Click += new System.EventHandler(this.AnalyzeButton_Click);
            // 
            // YIontxtBox
            // 
            this.YIontxtBox.Location = new System.Drawing.Point(24, 109);
            this.YIontxtBox.Name = "YIontxtBox";
            this.YIontxtBox.Size = new System.Drawing.Size(316, 22);
            this.YIontxtBox.TabIndex = 5;
            this.YIontxtBox.TextChanged += new System.EventHandler(this.YIontxtBox_TextChanged);
            // 
            // BIontxtBox
            // 
            this.BIontxtBox.Location = new System.Drawing.Point(24, 55);
            this.BIontxtBox.Name = "BIontxtBox";
            this.BIontxtBox.Size = new System.Drawing.Size(316, 22);
            this.BIontxtBox.TabIndex = 4;
            this.BIontxtBox.TextChanged += new System.EventHandler(this.BIontxtBox_TextChanged);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(21, 86);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(66, 17);
            this.label2.TabIndex = 3;
            this.label2.Text = "Y Ion File";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(21, 35);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(66, 17);
            this.label1.TabIndex = 2;
            this.label1.Text = "B Ion File";
            // 
            // YIonButton
            // 
            this.YIonButton.Location = new System.Drawing.Point(362, 107);
            this.YIonButton.Name = "YIonButton";
            this.YIonButton.Size = new System.Drawing.Size(69, 24);
            this.YIonButton.TabIndex = 1;
            this.YIonButton.Text = "Browse";
            this.YIonButton.UseVisualStyleBackColor = true;
            this.YIonButton.Click += new System.EventHandler(this.Button_2);
            // 
            // BIonButton
            // 
            this.BIonButton.Location = new System.Drawing.Point(362, 51);
            this.BIonButton.Name = "BIonButton";
            this.BIonButton.Size = new System.Drawing.Size(69, 24);
            this.BIonButton.TabIndex = 0;
            this.BIonButton.Text = "Browse";
            this.BIonButton.UseVisualStyleBackColor = true;
            this.BIonButton.Click += new System.EventHandler(this.Button_1);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1333, 490);
            this.Controls.Add(this.splitContainer1);
            this.Name = "Form1";
            this.Text = "Release 0.01";
            this.Load += new System.EventHandler(this.Form1_Load);
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).EndInit();
            this.splitContainer1.Panel1.ResumeLayout(false);
            this.splitContainer1.Panel1.PerformLayout();
            this.splitContainer1.Panel2.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer1)).EndInit();
            this.splitContainer1.ResumeLayout(false);
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataGridView dataGridView1;
        private System.Windows.Forms.Button BIonButton;
        private System.Windows.Forms.SplitContainer splitContainer1;
        private System.Windows.Forms.Button AnalyzeButton;
        private System.Windows.Forms.TextBox YIontxtBox;
        private System.Windows.Forms.TextBox BIontxtBox;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Button YIonButton;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Button FASTAButton;
        private System.Windows.Forms.TextBox FASTAtxtBox;
    }
}

