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

    gStyle->SetOptStat(1000010);

    TCanvas *cKanwa = new TCanvas("cKanwa","Canvas for plotting",1350,800);
    cKanwa->Divide(4,2);
    cKanwa->cd();


    cKanwa->cd(1);
    gPad->SetLogz();
    gPad->SetLogy();
    hst_xy->Draw("colz");

    cKanwa->cd(2);
    gPad->SetLogz();
    ptpz->Draw("colz");

    cKanwa->cd(3);
    gPad->SetLogz();
    ptpzg->Draw("colz");

    cKanwa->cd(4);
    gPad->SetLogz();
    pp->Draw("colz");

    cKanwa->cd(5);
    gPad->SetLogy();
    ph_angles->Draw();
    
    cKanwa->cd(6);
    gPad->SetLogz();
    ph_erg->Draw("colz");

    cKanwa->cd(7);
    gPad->SetLogy();

    kk->Draw();

    cKanwa->cd(8);
    w_int->Draw();
    gPad->SetLogy();
}