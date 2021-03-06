    void plotHistos(){
        
    TFile *f = new TFile("histos.root","read");

    TH1D *hst_xy = (TH1D*)f->Get("hst_xy");
    TH1D *ptpz = (TH1D*)f->Get("ptpz");
    TH1D *ptpzg = (TH1D*)f->Get("ptpzg");
    TH1D *pp = (TH1D*)f->Get("pp");
    TH1D *ph_angles = (TH1D*)f->Get("ph_angles");
    TH1D *ph_erg = (TH1D*)f->Get("ph_erg");
    TH1D *kk = (TH1D*)f->Get("kk");
    TH1D *w_int = (TH1D*)f->Get("w_int");

    gStyle->SetOptStat("");//1000010);

    TCanvas *plots = new TCanvas("plots","Moller/Bhabha Plots",1350,800);
    // TCanvas *moller = new TCanvas("moller","Photon Phase Space",1350,800);


    plots->Divide(4,2);
    plots->cd();


    plots->cd(1);
    gPad->SetLogz();
    gPad->SetLogy();
    hst_xy->Draw("colz");

    plots->cd(2);
    gPad->SetLogz();
    ptpz->Draw("colz");

    plots->cd(3);
    gPad->SetLogz();
    ptpzg->Draw("colz");

    plots->cd(4);
    gPad->SetLogz();
    pp->Draw("colz");

    plots->cd(5);
    gPad->SetLogy();
    ph_angles->Draw();
    
    plots->cd(6);
    gPad->SetLogz();
    ph_erg->Draw("colz");

    plots->cd(7);
    gPad->SetLogy();

    kk->Draw();

    plots->cd(8);
    w_int->Draw();
    gPad->SetLogy();
    
    // moller->cd();
    // gPad->SetLogz();
    // hst_xy->Draw("colz");

    FILE *oFile;
    oFile = fopen("output_ptpz.txt","w");
    for(int i = 0;i<ptpz->GetNbinsX();i++){
        for(int j = 0;j<ptpz->GetNbinsY();j++){
            // if(ptpz->GetBinContent(i,j)>0){
            fprintf(oFile,"%.3f\t%.3f\t%.3f\t%.3f\n",ptpz->GetXaxis()->GetBinCenter(i),
                ptpz->GetYaxis()->GetBinCenter(j),ptpz->GetBinContent(i,j),ptpz->GetBinError(i,j));
        // }
}
}
    fclose(oFile);

}
