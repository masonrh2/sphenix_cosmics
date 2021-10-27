//! Example program demonstrating parallel fitting with ROOT
//! by Andreas Zoglauer (andreas@megalibtoolkit.com)

//! Usage:
//! root
//! .L ParallelFitting.C++
//! ParallelFitting(4)
//! --> Like in all other multi-threaded examples, you have to quit root for a seconds run!

// Standard
#include <iostream>
#include <vector>
using namespace std;

// ROOT
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom.h"
#include "TThread.h"
#include "TMutex.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TStopwatch.h"

#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "HFitInterface.h"
#include "Fit/Fitter.h"


/******************************************************************************/

class ParallelFitter
{
public:
  //! Default constructor
  ParallelFitter() { m_NThreads = 4; }
  //! Default destructor
  ~ParallelFitter() {};
  
  //! Set the number of threads
  void SetThreads(unsigned int NThreads) { m_NThreads = NThreads; }
  
  //! Analyze what eveer needs to be analyzed...
  bool Analyze();
 
  //! In this function the parallel calibration happens
  bool CalibrateParallel(unsigned int ThreadID, bool IsMultiThreaded = false);


private:
  //! Input: some histograms
  vector<TH1D*> m_Histograms;
  //! Output: Mean of the fitted Gaussian
  vector<double> m_Means;
  
  //! Number of threads
  unsigned int m_NThreads;
  //! Storing the threads:
  vector<TThread*> m_Threads;
  //! Storing a flag that the thread is running
  vector<bool> m_ThreadIsInitialized;
  //! Storing a flag that the thread is finished
  vector<bool> m_ThreadIsFinished;
  //! ID of the next item to be processed
  unsigned int m_ThreadNextItem;

};


/******************************************************************************
 * This function encapsulates the information transfered into the thread
 */
class ThreadCaller
{
 public:
  //! Standard constructor
  ThreadCaller(ParallelFitter* M, unsigned int ThreadID) {
    m_Fitter = M;
    m_ThreadID = ThreadID;
  }

  //! Return the calling class
  ParallelFitter* GetThreadCaller() { return m_Fitter; }
  //! Return the thread ID
  unsigned int GetThreadID() { return m_ThreadID; }

 private:
  //! Store the calling class for retrieval
  ParallelFitter* m_Fitter;
  //! ID of the worker thread
  unsigned int m_ThreadID;
};


/******************************************************************************
 * Thread entry point for the parallel calibration
 */
void CallParallelCalibrationThread(void* Address)
{
  ParallelFitter* P = ((ThreadCaller*) Address)->GetThreadCaller();
  P->CalibrateParallel(((ThreadCaller*) Address)->GetThreadID(), true);
}


/******************************************************************************
 * The fit function - Gaus
 */
class Gaus : public ROOT::Math::IParamFunction { 
public:
  void SetParameters(const double* p) { std::copy(p,p+NPar(),fp);}
  const double* Parameters() const { return fp; }
  ROOT::Math::IGenFunction* Clone() const { 
    Gaus* g = new Gaus(); 
    g->SetParameters(fp);
    return g;
  };
  unsigned int NPar() const { return 4; }
private:
  double DoEvalPar(double x, const double* par) const { 
    double fitval = par[0];
    double arg = 0;

    if (par[3] != 0) arg = (x - par[2])/par[3];
    fitval += par[1]*TMath::Exp(-0.5*arg*arg);

    return fitval;
  }
  double fp[4];  
};


/******************************************************************************
 * Do whatever analysis is necessary
 */
bool ParallelFitter::Analyze()
{
  // Sanity checks:
  // We will be using Minuit2 for fitting since it is thread safe, so chekc it is there
  if (gROOT->GetPluginManager()->FindHandler("ROOT::Math::Minimizer", "Minuit2") == 0)  {
    cout<<"We need Minuit2 for parallel fitting!"<<endl;
    return false;
  }
  
  
  // Create the list of histograms:
  unsigned int Histograms = 5000;
  unsigned int Events = 5000;
  cout<<"Creating "<<Histograms<<" histograms with "<<Events<<" events for parallel fitting..."<<endl;
  for (unsigned int h = 0; h < Histograms; ++h) {
    TH1D* Hist = new TH1D("", "Hist", 100, -10, 10);
    for (unsigned int i = 0; i < Events; ++i) {
      Hist->Fill(gRandom->Gaus(0, 1));
    }
    m_Histograms.push_back(Hist);
    m_Means.push_back(-999);
  }
    
  m_ThreadNextItem = 0;
  m_Threads.resize(m_NThreads);
  m_ThreadIsInitialized.resize(m_NThreads);
  m_ThreadIsFinished.resize(m_NThreads);
    
  TStopwatch StopWatch;
  
  if (m_NThreads > 1) {
    
    // Start threads
    for (unsigned int t = 0; t < m_NThreads; ++t) {
      TString Name = "Calibration thread #";
      Name += t;
      cout<<"Creating thread: "<<Name<<endl;
      TThread* Thread = new TThread(Name, (void(*) (void *)) &CallParallelCalibrationThread, (void*) new ThreadCaller(this, t));
      m_Threads[t] = Thread;
      m_ThreadIsInitialized[t] = false;
      m_ThreadIsFinished[t] = false;

      Thread->Run();
      
      // Wait until thread is initialized:
      while (m_ThreadIsInitialized[t] == false && m_ThreadIsFinished[t] == false) {
        // Sleep for a while...
        TThread::Sleep(0, 10000000);
      }    

      cout<<Name<<" is running"<<endl;
    }
    
    bool ThreadsAreRunning = true;
    while (ThreadsAreRunning == true) {

      // Sleep for a while...
      TThread::Sleep(0, 10000000);
      
      ThreadsAreRunning = false;
      for (unsigned int t = 0; t < m_NThreads; ++t) {
        if (m_ThreadIsFinished[t] == false) {
          ThreadsAreRunning = true;
          break;
        }
      }
    }

    // None of the threads are running any more --- kill them
    for (unsigned int t = 0; t < m_NThreads; ++t) {
      m_Threads[t]->Kill();
      m_Threads[t] = 0;
    }
    
    cout<<"All threads have finished"<<endl;
  }
  // Non-threaded mode
  else {
    CalibrateParallel(0);
  }
  
  cout<<"Elapsed time in parallel section: "<<StopWatch.RealTime()<<endl;
  
  // Create a histogram of all means
  TH1D* Means = new TH1D("", "Distribution of all fitted Gaussian means", 100, -0.1, 0.1);
  for (unsigned int i = 0; i < m_Means.size(); ++i) {
    Means->Fill(m_Means[i]);  
  }
  
  TCanvas* C = new TCanvas();
  C->cd();
  Means->Draw();
  C->Update();
  
  // Delete the histograms
  for (unsigned int i = 0; i < m_Histograms.size(); ++i) {
    delete m_Histograms[i]; 
  }
  m_Histograms.clear();
  m_Means.clear();
  
  return true;
}


/******************************************************************************
 * Perform the parallel calibration 
 */
bool ParallelFitter::CalibrateParallel(unsigned int ThreadID, bool IsMultiThreaded)
{
  if (IsMultiThreaded == true) {
    cout<<"Calibrate parallel for thread #"<<ThreadID<<" called"<<endl;
  }
  
  m_ThreadIsInitialized[ThreadID] = true;
  
  Gaus FitFunction;
  ROOT::Fit::Fitter TheFitter; 
  TheFitter.Config().SetMinimizer("Minuit2");
  TheFitter.SetFunction(FitFunction);

  while (true) {
    TThread::Lock();
    unsigned int ID = m_ThreadNextItem;
    ++m_ThreadNextItem;
    TThread::UnLock();
    
    if (ID >= m_Histograms.size()) break;
    
    //cout<<"Thread #"<<ThreadID<<": Taking care of histogram #"<<ID<<endl;
    
    ROOT::Fit::BinData d; 
    ROOT::Fit::FillData(d, m_Histograms[ID]);
    
    TheFitter.Config().ParSettings(0).SetName("Offset");
    TheFitter.Config().ParSettings(0).SetValue(0);
    TheFitter.Config().ParSettings(0).SetLimits(0, 0);
    TheFitter.Config().ParSettings(1).SetName("Height");
    TheFitter.Config().ParSettings(1).SetValue(1);
    TheFitter.Config().ParSettings(1).SetLimits(0, 10000);
    TheFitter.Config().ParSettings(2).SetName("Mean");
    TheFitter.Config().ParSettings(2).SetValue(1);
    TheFitter.Config().ParSettings(2).SetLimits(-10, 10);
    TheFitter.Config().ParSettings(3).SetName("Sigma");
    TheFitter.Config().ParSettings(3).SetValue(1);
    TheFitter.Config().ParSettings(3).SetLimits(0.1, 10);

    bool ReturnCode = false;
    ReturnCode = TheFitter.Fit(d);
   
    if (ReturnCode == true) {
      m_Means[ID] = TheFitter.Result().Parameter(2);
    }
  }
    
  m_ThreadIsFinished[ThreadID] = true;

  if (IsMultiThreaded == true) {
    cout<<"Thread #"<<ThreadID<<" has finished processing"<<endl;
  }
  
  return true;
}


/******************************************************************************
 * Main program
 */
int test_mt(int Threads = 4)
//int main()
{
  //int Threads = 4;
  
  ParallelFitter P;
  P.SetThreads(Threads);
  P.Analyze();

  cout<<"Program finished!"<<endl;

  return 0;
}