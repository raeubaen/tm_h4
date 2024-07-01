#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <TRandom.h>
#include <TH2F.h>

template <int dim>
struct Vector
{
    std::array<double, dim> components; 
};


using Vector2D = Vector<2>;
using Vector3D = Vector<3>;


// subtract
template <int dim>
Vector<dim> operator-(const Vector<dim> &u, const Vector<dim> &v) {
    Vector<dim> result(u);
    for (int i = 0; i < dim; ++i)
        result.components[i] -= v.components[i];
    return result;
}

// add
template <int dim>
Vector<dim> operator+(const Vector<dim> &u, const Vector<dim> &v) {
    Vector<dim> result(u);
    for (int i = 0; i < dim; ++i)
        result.components[i] += v.components[i];
    return result;
}

// negate
template <int dim>
Vector<dim> operator-(const Vector<dim> &u) {
    Vector<dim> result;
    for (int i = 0; i < dim; ++i)
        result.components[i] = -u.components[i];
    return result;
}

// scalar product
template <int dim>
double operator*(const Vector<dim> &u, const Vector<dim> &v) {
    double result = 0;
    for (int i = 0; i < dim; ++i)
        result += u.components[i] * v.components[i];
    return result;
}

// scale
template <int dim>
Vector<dim> operator*(const Vector<dim> &u, const double s) {
    Vector<dim> result(u);
    for (int i = 0; i < dim; ++i)
        result.components[i] *= s;
    return result;
}

// scale
template <int dim>
Vector<dim> operator*(const double s, const Vector<dim> &u) {
    return u*s;
}


// ostream
template <int dim>
std::ostream& operator<< (std::ostream& out, const Vector<dim> &u) {
    out << "(";
    for (auto c : u.components)
        out << std::setw(15) << c ;
    out << ")";
    return out;
}


std::pair<Vector3D, Vector3D>
shortest_connection_segment_to_segment(const Vector3D A, const Vector3D B, const Vector3D C, const Vector3D D)
{
    Vector3D u = B - A;
    Vector3D v = D - C;
    Vector3D w = A - C;

    double    a = u*u;         // always >= 0
    double    b = u*v;
    double    c = v*v;         // always >= 0
    double    d = u*w;
    double    e = v*w;
    double    sc, sN, sD = a*c - b*b;  // sc = sN / sD, default sD = D >= 0
    double    tc, tN, tD = a*c - b*b;  // tc = tN / tD, default tD = D >= 0
    double    tol = 1e-15;
    // compute the line parameters of the two closest points
    if (sD < tol) {            // the lines are almost parallel
        sN = 0.0;              // force using point A on segment AB
        sD = 1.0;              // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                     // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0;          // compute shortest connection of A to segment CD
            tN = e;
            tD = c;
        }
        else if (sN > sD) {    // sc > 1  => the s=1 edge is visible
            sN = sD;           // compute shortest connection of B to segment CD
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;             
        // recompute sc for this edge
        if (-d < 0.0)          // compute shortest connection of C to segment AB
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)  // compute shortest connection of D to segment AB
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d +  b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (fabs(sN) < tol ? 0.0 : sN / sD);
    tc = (fabs(tN) < tol ? 0.0 : tN / tD);

    Vector3D P1 = A + (sc * u);
    Vector3D P2 = C + (tc * v);

    return {P1, P2};   // return the closest distance
}


void toyVertex()
{
   // Get old file, old tree and set top branch address
   TString filename = "prova_reconuova.root";

   TFile oldfile(filename);
   TTree *oldtree;
   oldfile.GetObject("tree", oldtree);

   const auto nentries = oldtree->GetEntries();

   TFile newfile("prova_reconuova_2.root", "recreate");
   auto newtree = oldtree->CloneTree(0);
   double vze[1], vzp[1];
   newtree->Branch("vze", vze, "vze/D");
   newtree->Branch("vzp", vzp, "vzp/D");
   float reco_x1_e, reco_theta_e, reco_phi_e, reco_y1_e, reco_x2_e, reco_y2_e, reco_x1_p, reco_y1_p, reco_x2_p, reco_y2_p, reco_theta_p, reco_phi_p;
   oldtree->SetBranchAddress("reco_x1_f_e", &reco_x1_e);
   oldtree->SetBranchAddress("reco_x1_s_e", &reco_x2_e);
   oldtree->SetBranchAddress("reco_y1_f_e", &reco_y1_e);
   oldtree->SetBranchAddress("reco_y1_s_e", &reco_y2_e);
   oldtree->SetBranchAddress("reco_x1_f_p", &reco_x1_p);
   oldtree->SetBranchAddress("reco_x1_s_p", &reco_x2_p);
   oldtree->SetBranchAddress("reco_y1_f_p", &reco_y1_p);
   oldtree->SetBranchAddress("reco_y1_s_p", &reco_y2_p);
   oldtree->SetBranchAddress("reco_theta_e", &reco_theta_e);
   oldtree->SetBranchAddress("reco_phi_e", &reco_phi_e);
   oldtree->SetBranchAddress("reco_theta_p", &reco_theta_p);
   oldtree->SetBranchAddress("reco_phi_p", &reco_phi_p);


   for (int i=0; i<int(nentries); i++) {
      if (i%1000==0) cout << i << endl;
      oldtree->GetEntry(i);

      Vector3D e1 = {reco_x1_e-2000*sin(reco_theta_e)*cos(reco_phi_e), reco_y1_e-2000*sin(reco_theta_e)*sin(reco_phi_e), -2000*cos(reco_theta_e)};
      Vector3D e2 = {reco_x1_e, reco_y1_e, 0};
      Vector3D p1 = {reco_x1_p-2000*sin(reco_theta_p)*cos(reco_phi_p), reco_y1_p-2000*sin(reco_theta_p)*sin(reco_phi_p), -2000*cos(reco_theta_p)};
      Vector3D p2 = {reco_x1_p, reco_y1_p, 0};

      auto [VE, VP] = shortest_connection_segment_to_segment(e1, e2, p1, p2);
      vze[0] = VE.components[2];
      vzp[0] = VP.components[2];
      newtree->Fill();
   }

   newtree->Print();
   newfile.Write();
}
