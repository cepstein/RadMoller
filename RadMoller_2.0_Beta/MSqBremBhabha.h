double RadMoller_Gen::Mh2b(TLorentzVector *p1, TLorentzVector *p2,\
 TLorentzVector *q1, TLorentzVector *q2, TLorentzVector *k){
  kp1 = k->Dot(*p1);
  kp2 = k->Dot(*p2);
  kp3 = k->Dot(*q1);
  p1p2 = p1->Dot(*p2);
  p1p3 = p1->Dot(*q1);

  return -(pow(ec,6)*(((512*(3*pow(me,6) + 
               pow(me,4)*(-kp1 - kp2 + 2*(kp1 + kp2 - kp3) + 4*p1p2) + 
               kp1*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
               (kp1 + kp2 - 2*kp3)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*p1p3 + 
               kp2*(p1p3*p1p3) + 
               p1p2*((kp2 - kp3)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  (kp1 - kp3)*p1p3) + 
               me*me*((kp1 + kp2 - 2*kp3)*p1p2 + p1p2*p1p2 + 
                  (kp1 + 2*kp2 - kp3)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  (2*kp1 + kp2 - kp3 - 
                     2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*p1p3)))/
           pow(me*me + p1p2,2) + 
          (128*(30*pow(me,6) + 
               16*((kp2 - kp3)*(p1p2 - p1p3) + 
                  kp1*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
               pow(me,4)*(-4*kp1 - kp2 + 4*(kp1 + kp2 - kp3) + 
                  7*(me*me) + 7*p1p2 + 
                  2*(-2*kp1 + 5*kp2 - 10*kp3 - me*me + 18*p1p2 + 
                     10*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     16*p1p3) - 10*(-kp1 + me*me + p1p2 - p1p3) - 
                  4*(kp1 - kp3 + p1p3)) + 
               me*me*(-2*kp1*(me*me) + 3*kp2*(me*me) + 4*kp3*(me*me) - 
                  pow(me,4) + 8*(p1p2*p1p2) + 
                  36*kp1*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  36*kp2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  36*kp3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  8*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  4*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                  2*(8*kp2 - 8*kp3 + 3*(me*me) + 
                     6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*p1p3 + 
                  8*(p1p3*p1p3) - 
                  p1p2*(-16*kp2 + 16*kp3 + me*me - 
                     12*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     16*p1p3))))/((me*me + p1p2)*(2*(me*me) - 2*p1p3)) + 
          (128*(56*pow(me,6) + 
               me*me*(-6*kp1*(me*me) + kp2*(me*me) + 12*kp3*(me*me) + 
                  2*pow(me,4) + 
                  p1p2*(-32*kp1 - 16*kp2 + 16*kp3 + me*me + 
                     24*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                  24*kp1*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  24*kp2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  40*kp3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  4*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  8*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                  2*(8*kp1 + 16*kp2 - 8*kp3 + 13*(me*me) - 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*p1p3 + 
                  16*(p1p3*p1p3)) + 
               pow(me,4)*(18*(me*me) + 24*p1p2 + 
                  12*(-kp1 + me*me + p1p2 - p1p3) - 
                  49*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  60*p1p3 - 3*
                   (kp1 - kp3 + me*me + 2*(-kp1 + me*me + p1p2 - p1p3) - 
                     2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + p1p3)) \
+ 16*(-(kp3*(p1p2*p1p2)) + p1p2*
                   ((kp1 + 2*kp2 - kp3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     (kp1 + kp2)*p1p3) + 
                  (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (kp1*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     (-kp2 + kp3)*p1p3))))/pow(2*(me*me) - 2*p1p3,2))/
        (4.*pow(kp1 + kp2 - kp3,2)) - 
       (32*(-192*pow(me,6) + 31*(p1p2*p1p2)*(-kp1 - kp2 + p1p2) + 
            pow(me,4)*(27*kp1 - 39*(kp1 + kp2 - kp3) - 26*kp3 + 
               181*(me*me) - 6*p1p2 + 56*(-kp1 - kp2 + p1p2) + 
               57*(-kp1 + me*me + p1p2 - p1p3) + 
               32*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 80*p1p3 + 
               67*(kp1 - kp3 + p1p3)) + 
            p1p2*((-kp1 + me*me + p1p2 - p1p3)*
                (16*(me*me) - 63*(-kp1 - kp2 + p1p2) + 
                  66*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 16*p1p3) \
+ (-kp1 - kp2 + p1p2)*(104*(me*me) + 31*(-kp1 - kp2 + p1p2) - 
                  31*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  31*(kp1 - kp3 + p1p3)) + 
               p1p3*(-47*(-kp1 - kp2 + p1p2) + 6*(kp1 - kp3 + p1p3))) + 
            2*(5*(me*me)*(-kp1 - kp2 + p1p2)*(-kp2 + me*me + p1p2) - 
               33*pow(-kp1 + me*me + p1p2 - p1p3,2)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
               3*(p1p3*p1p3)*(kp1 - kp3 + p1p3) + 
               (-kp1 + me*me + p1p2 - p1p3)*
                (-29*(me*me)*(-kp1 - kp2 + p1p2) - 
                  33*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                  8*(me*me - 2*(-kp1 - kp2 + p1p2))*(kp1 - kp3 + p1p3) + 
                  (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (116*(me*me) + 49*(-kp1 - kp2 + p1p2) - 
                     33*(kp1 - kp3 + p1p3))) + 
               p1p3*((-kp1 - kp2 + p1p2)*
                   (-13*(me*me) + 
                     8*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                  (12*(me*me) + 11*(-kp1 - kp2 + p1p2) - 
                     3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                   (kp1 - kp3 + p1p3) - 3*pow(kp1 - kp3 + p1p3,2) + 
                  (-kp1 + me*me + p1p2 - p1p3)*
                   (16*(me*me) - 
                     41*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     11*(kp1 - kp3 + p1p3)))) + 
            me*me*(117*pow(me,4) + p1p2*p1p2 + 
               156*(me*me)*(-kp1 - kp2 + p1p2) + 
               61*pow(-kp1 - kp2 + p1p2,2) - 
               154*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
               97*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) + 
               4*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
               77*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
               59*(-kp1 - kp2 + p1p2)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
               2*(-kp1 + me*me + p1p2 - p1p3)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
               18*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
               42*(p1p3*p1p3) - 95*(me*me)*(kp1 - kp3 + p1p3) - 
               57*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
               38*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                (kp1 - kp3 + p1p3) - 20*pow(kp1 - kp3 + p1p3,2) + 
               p1p2*(88*(me*me) + 78*(-kp1 - kp2 + p1p2) - 
                  21*(-kp1 + me*me + p1p2 - p1p3) + 
                  17*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 57*p1p3 + 
                  19*(kp1 - kp3 + p1p3)) - 
               p1p3*(18*(me*me) + 35*(-kp1 - kp2 + p1p2) + 
                  22*(-kp1 + me*me + p1p2 - p1p3) + 
                  76*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  78*(kp1 - kp3 + p1p3)) - 
               2*(p1p2*(21*(-kp1 - kp2 + p1p2) + 
                     8*(-kp1 + me*me + p1p2 - p1p3)) - 
                  13*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) + 
                  5*(-kp1 - kp2 + p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  50*(-kp1 + me*me + p1p2 - p1p3)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  5*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                  8*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                  p1p3*(-5*(-kp1 - kp2 + p1p2) + 
                     8*(-kp1 + me*me + p1p2 - p1p3) + 6*(kp1 - kp3 + p1p3))\
))))/(kp2*kp2*pow(2*(me*me) - 2*p1p3,2)) + 
       ((256*(6*pow(me,6) + 
               pow(me,4)*(-kp1 + kp2 + 
                  2*(kp2 + 4*(-kp1 - kp2 + p1p2))) + 
               2*(p1p2*p1p2)*p1p3 + 
               (-kp1 + kp3 + 2*(me*me) - p1p3)*(p1p3*p1p3) - 
               (-kp1 + me*me + p1p2 - p1p3)*
                (-(me*me*(-kp1 - kp2 + p1p2)) + 
                  (-kp1 - kp2 - kp3 - me*me + 2*p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                  (me*me + 2*(-kp1 - kp2 + p1p2))*(kp1 - kp3 + p1p3)) + 
               me*me*(p1p2*p1p2 + 4*(me*me)*(-kp1 - kp2 + p1p2) + 
                  3*pow(-kp1 - kp2 + p1p2,2) - 
                  p1p2*(-kp1 + me*me + p1p2) + 
                  3*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
                  (-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) + 
                  me*me*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  3*(-kp1 - kp2 + p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  4*(-kp1 + me*me + p1p2 - p1p3)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  2*(p1p3*p1p3) + me*me*(kp1 - kp3 + p1p3) - 
                  3*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                  4*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                  (-kp1 + me*me + p1p2 - p1p3)*
                   (-kp2 - kp3 + p1p2 + p1p3) - 
                  p1p3*(2*kp1 + 2*kp2 - kp3 - me*me - 2*p1p2 + 
                     2*(-kp1 + me*me + p1p2 - p1p3) + p1p3) + 
                  p1p2*(kp2 + me*me - p1p2 + 4*(-kp1 - kp2 + p1p2) + 
                     3*(-kp1 + me*me + p1p2 - p1p3) + 3*p1p3) + 
                  p1p3*(-kp1 - kp2 + 5*(me*me) + p1p2 - 
                     6*(-kp1 + me*me + p1p2 - p1p3) - 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     4*(kp1 - kp3 + p1p3))) + 
               p1p2*((-kp1 + me*me + p1p2 - p1p3)*
                   (-kp1 - kp2 + kp3 + 2*(me*me) + p1p2 + 
                     2*(-kp1 - kp2 + p1p2) - p1p3) + 
                  p1p3*(5*(me*me) + 2*(-kp1 - kp2 + p1p2) - 
                     4*(-kp1 + me*me + p1p2 - p1p3) - 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     3*(kp1 - kp3 + p1p3))) + 
               p1p3*(-(me*me*(-kp1 - kp2 + p1p2)) + 
                  2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                  (-kp1 - kp2 + 2*(me*me) + p1p2 - 
                     3*(-kp1 + me*me + p1p2 - p1p3))*(kp1 - kp3 + p1p3) \
+ pow(kp1 - kp3 + p1p3,2) + 
                  (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (-5*(me*me) - 2*(-kp1 - kp2 + p1p2) + 
                     3*(-kp1 + me*me + p1p2 - p1p3) + 
                     3*(kp1 - kp3 + p1p3)))))/(kp2*kp2) + 
          (256*(6*pow(me,6) + 
               pow(me,4)*(kp1 - kp2 + 
                  2*(kp1 + 4*(-kp1 - kp2 + p1p2))) + 
               2*(p1p2*p1p2)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
               me*me*(-kp1 - kp2 + p1p2)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
               2*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
               (-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
               pow(-kp1 + me*me + p1p2 - p1p3,2)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
               2*(me*me)*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                 2) - (-kp1 + me*me + p1p2 - p1p3)*
                pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
               (-kp1 + kp3 + 
                  2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - p1p3)*
                (p1p3*p1p3) + 
               me*me*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
               me*me*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) - 
               2*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)*
                (kp1 - kp3 + p1p3) + 
               3*(-kp1 + me*me + p1p2 - p1p3)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                (kp1 - kp3 + p1p3) + 
               p1p2*((kp1 - kp3 + p1p3)*
                   (me*me + 2*(-kp1 - kp2 + p1p2) + p1p3) + 
                  (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (5*(me*me) + 2*(-kp1 - kp2 + p1p2) - 
                     3*(-kp1 + me*me + p1p2 - p1p3) - 4*p1p3 - 
                     4*(kp1 - kp3 + p1p3))) - 
               p1p3*((-kp1 - kp2 - kp3 - me*me + 2*p1p2)*
                   (kp1 - kp3 + p1p3) + 
                  (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (5*(me*me) + 2*(-kp1 - kp2 + p1p2) - 
                     3*(-kp1 + me*me + p1p2 - p1p3) - 
                     3*(kp1 - kp3 + p1p3))) + 
               me*me*(p1p2*p1p2 + 4*(me*me)*(-kp1 - kp2 + p1p2) + 
                  3*pow(-kp1 - kp2 + p1p2,2) - 
                  p1p2*(-kp2 + me*me + p1p2) + 
                  me*me*(-kp1 + me*me + p1p2 - p1p3) - 
                  3*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) + 
                  5*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  2*(-kp1 - kp2 + p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  4*(-kp1 + me*me + p1p2 - p1p3)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                  (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*p1p3 + 
                  3*(me*me)*(kp1 - kp3 + p1p3) - 
                  3*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) - 
                  8*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (kp1 - kp3 + p1p3) + 
                  p1p3*(me*me - 3*(-kp1 - kp2 + p1p2) - 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     4*(kp1 - kp3 + p1p3)) + 
                  p1p2*(kp1 + me*me - p1p2 + 4*(-kp1 - kp2 + p1p2) + 
                     3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     3*(kp1 - kp3 + p1p3)))))/(kp1*kp1) - 
          (256*(4*pow(me,6) + 
               2*(p1p2*p1p2)*(-kp1 - kp2 + kp3 + me*me + p1p2 + 
                  2*(-kp1 - kp2 + p1p2)) - 
               me*me*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) + 
               me*me*(-kp1 - kp2 + p1p2)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
               6*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
               2*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
               2*pow(-kp1 + me*me + p1p2 - p1p3,2)*
                (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
               4*(-kp1 + me*me + p1p2 - p1p3)*
                pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
               me*me*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
               2*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                (kp1 - kp3 + p1p3) + 
               2*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                (kp1 - kp3 + p1p3) + 
               2*pow(me,4)*(-2*(me*me) + 8*p1p2 - 
                  2*(-kp1 - kp2 + p1p2) + 
                  2*(-kp1 + me*me + p1p2 - p1p3) + 2*(kp1 - kp3 + p1p3)) \
+ 2*(p1p3*p1p3)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3 + 
                  2*(kp1 - kp3 + p1p3)) + 
               p1p3*(me*me*(-kp1 - kp2 + p1p2) + 
                  2*(me*me)*(-kp1 + me*me + p1p2 - p1p3) + 
                  2*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                  6*(me*me)*(kp1 - kp3 + p1p3) - 
                  2*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                  2*pow(kp1 - kp3 + p1p3,2)) + 
               p1p2*(8*(me*me)*(-kp1 - kp2 + p1p2) + 
                  4*pow(-kp1 - kp2 + p1p2,2) + 
                  (me*me - 2*(-kp1 - kp2 + p1p2) - 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                   (-kp1 + me*me + p1p2 - p1p3) + 
                  5*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  2*(-kp1 - kp2 + p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                  2*(p1p3*p1p3) + me*me*(kp1 - kp3 + p1p3) - 
                  2*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                  2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (kp1 - kp3 + p1p3) + 
                  p1p3*(5*(me*me) - 2*(-kp1 - kp2 + p1p2) - 
                     2*(-kp1 + me*me + p1p2 - p1p3) - 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     4*(kp1 - kp3 + p1p3))) + 
               me*me*(-(p1p2*(-kp1 - kp2 + 2*(me*me) + 2*p1p2)) - 
                  (-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) - 
                  (-kp1 - kp2 + p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  (-kp1 - kp2 + p1p2 + 2*(-kp1 + me*me + p1p2 - p1p3) + 
                     2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*p1p3 - 
                  (-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                  2*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) - 
                  2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (kp1 - kp3 + p1p3) + 
                  2*(p1p2*p1p2 - pow(-kp1 - kp2 + p1p2,2) + 
                     5*(me*me)*(-kp1 + me*me + p1p2 - p1p3) + 
                     5*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) - 
                     3*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
                     4*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     p1p3*p1p3 - me*me*(kp1 - kp3 + p1p3) + 
                     2*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                     (-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                     p1p2*(kp1 - kp3 - 2*(me*me) + 
                        4*(-kp1 - kp2 + p1p2) + 
                        4*(-kp1 + me*me + p1p2 - p1p3) + 2*p1p3) + 
                     p1p3*(-kp3 + me*me - 
                        4*(-kp1 + me*me + p1p2 - p1p3) + p1p3 - 
                        2*(kp1 - kp3 + p1p3))))))/(kp1*kp2))/
        (4.*pow(-kp1 - kp2 + me*me + p1p2,2)) + 
       128*((42*pow(me,6) + 
             pow(me,4)*(-2*(me*me) - 3*p1p2 - 10*(-kp1 - kp2 + p1p2) - 
                16*(-kp1 + me*me + p1p2 - p1p3) + 
                22*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 26*p1p3 - 
                56*(kp1 - kp3 + p1p3)) + 
             4*(p1p2*p1p2*(-kp1 - kp2 + p1p2 + 
                   4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                me*me*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) + 
                me*me*(-kp1 - kp2 + p1p2)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                4*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                6*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                2*pow(-kp1 + me*me + p1p2 - p1p3,2)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                4*(me*me)*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                  2) - 2*(-kp1 + me*me + p1p2 - p1p3)*
                 pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*(p1p3*p1p3) + 
                3*(me*me)*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                4*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)*
                 (kp1 - kp3 + p1p3) + 
                2*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                 (kp1 - kp3 + p1p3) + 
                2*(-kp1 + me*me + p1p2 - p1p3)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                 (kp1 - kp3 + p1p3) + 
                p1p2*((-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (10*(me*me) + 7*(-kp1 - kp2 + p1p2) - 
                      6*(-kp1 + me*me + p1p2 - p1p3) - 8*p1p3 - 
                      4*(kp1 - kp3 + p1p3)) + 
                   (-kp1 - kp2 + p1p2)*
                    (-kp2 + me*me + 3*(kp1 - kp3 + p1p3))) + 
                p1p3*((-kp1 - kp2 + p1p2)*
                    (me*me - 4*(kp1 - kp3 + p1p3)) + 
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (-10*(me*me) - 8*(-kp1 - kp2 + p1p2) + 
                      6*(-kp1 + me*me + p1p2 - p1p3) + 
                      4*(kp1 - kp3 + p1p3)))) + 
             me*me*(-8*(p1p2*p1p2) - 
                p1p2*(11*(me*me) + 32*(-kp1 - kp2 + p1p2) - 
                   16*(-kp1 + me*me + p1p2 - p1p3) - 
                   32*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   28*p1p3 + 4*(kp1 - kp3 + p1p3)) - 
                4*((-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) + 
                   2*p1p2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   3*(-kp1 - kp2 + p1p2)*
                    (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                   (-kp1 - kp2 + p1p2 - 
                      2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*p1p3 - 
                   (-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                   2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (kp1 - kp3 + p1p3)) - 
                2*(2*pow(me,4) + 10*(me*me)*(-kp1 - kp2 + p1p2) - 
                   4*pow(-kp1 - kp2 + p1p2,2) + 
                   4*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
                   12*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   16*(-kp1 - kp2 + p1p2)*
                    (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                   4*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                   10*(p1p3*p1p3) + me*me*(kp1 - kp3 + p1p3) + 
                   2*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                   6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (kp1 - kp3 + p1p3) - 6*pow(kp1 - kp3 + p1p3,2) - 
                   (-kp1 + me*me + p1p2 - p1p3)*
                    (15*(me*me) + 16*(-kp1 - kp2 + p1p2) - 
                      16*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      2*(kp1 - kp3 + p1p3)) + 
                   p1p3*(-6*(me*me) - 22*(-kp1 - kp2 + p1p2) + 
                      14*(-kp1 + me*me + p1p2 - p1p3) + 
                      22*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      4*(kp1 - kp3 + p1p3)))))/
           (4.*(kp1*kp1)*pow(2*(me*me) - 2*(kp1 - kp3 + p1p3),2)) + 
          (14*pow(me,6) + pow(me,4)*
              (-86*kp1 + 14*(kp1 + kp2 - kp3) + 67*kp3 + 80*(me*me) + 
                26*p1p2 + 30*(-kp1 - kp2 + p1p2) - 
                82*(-kp1 + me*me + p1p2 - p1p3) - 
                40*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 78*p1p3 - 
                28*(kp1 - kp3 + p1p3)) + 
             me*me*(53*pow(me,4) + 26*(p1p2*p1p2) + 
                25*(me*me)*(-kp1 - kp2 + p1p2) + 
                36*pow(-kp1 - kp2 + p1p2,2) + 
                24*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
                4*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
                205*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                86*(-kp1 - kp2 + p1p2)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                46*(-kp1 + me*me + p1p2 - p1p3)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                114*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                72*(p1p3*p1p3) - 28*(me*me)*(kp1 - kp3 + p1p3) - 
                18*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                6*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                112*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                 (kp1 - kp3 + p1p3) + 14*pow(kp1 - kp3 + p1p3,2) + 
                2*p1p2*(23*(me*me) + 7*(-kp1 - kp2 + p1p2) - 
                   11*(-kp1 + me*me + p1p2 - p1p3) - 
                   62*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                   31*p1p3 - 12*(kp1 - kp3 + p1p3)) + 
                p1p3*(55*(me*me) + 68*(-kp1 - kp2 + p1p2) - 
                   92*(-kp1 + me*me + p1p2 - p1p3) + 
                   42*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   10*(kp1 - kp3 + p1p3)) - 
                2*(p1p2*(14*(-kp1 - kp2 + p1p2) - 
                      4*(-kp1 + me*me + p1p2 - p1p3) + 
                      8*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                   3*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) - 
                   23*(-kp1 - kp2 + p1p2)*
                    (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   28*(-kp1 + me*me + p1p2 - p1p3)*
                    (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   7*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                   4*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                   8*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (kp1 - kp3 + p1p3) + 
                   p1p3*(7*(-kp1 - kp2 + p1p2) - 
                      4*(-kp1 + me*me + p1p2 - p1p3) - 
                      16*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      4*(kp1 - kp3 + p1p3)))) + 
             2*(2*(p1p2*p1p2)*
                 (-3*(-kp1 - kp2 + p1p2) + 
                   8*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                19*(me*me)*(-kp1 - kp2 + p1p2)*
                 (-kp1 + me*me + p1p2 - p1p3) - 
                15*(me*me)*(-kp1 - kp2 + p1p2)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                72*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                42*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                26*pow(-kp1 + me*me + p1p2 - p1p3,2)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                18*(-kp1 + me*me + p1p2 - p1p3)*
                 pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                me*me*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                4*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                 (kp1 - kp3 + p1p3) - 
                16*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)*
                 (kp1 - kp3 + p1p3) + 
                8*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                 (kp1 - kp3 + p1p3) + 
                26*(-kp1 + me*me + p1p2 - p1p3)*
                 (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                 (kp1 - kp3 + p1p3) + 
                2*(p1p3*p1p3)*
                 (12*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   5*(kp1 - kp3 + p1p3)) + 
                p1p3*(15*(me*me)*(-kp1 - kp2 + p1p2) + 
                   16*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                   24*(me*me)*(kp1 - kp3 + p1p3) - 
                   6*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                   10*pow(kp1 - kp3 + p1p3,2) - 
                   2*(-kp1 + me*me + p1p2 - p1p3)*
                    (kp1 - kp3 + 2*(me*me) - 
                      29*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + p1p3) \
- 2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (32*(me*me) + 20*(-kp1 - kp2 + p1p2) - 
                      7*(kp1 - kp3 + p1p3))) - 
                2*p1p2*(-9*(me*me)*(-kp1 - kp2 + p1p2) - 
                   5*pow(-kp1 - kp2 + p1p2,2) + 
                   (2*(me*me) + 5*(-kp1 - kp2 + p1p2) + 
                      21*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                    (-kp1 + me*me + p1p2 - p1p3) + 
                   8*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                   (-kp1 - kp2 + p1p2)*p1p3 + 
                   (-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                   5*p1p3*(kp1 - kp3 + p1p3) + 
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (-20*(me*me) - 7*(-kp1 - kp2 + p1p2) + 20*p1p3 + 
                      8*(kp1 - kp3 + p1p3)))))/
           (4.*kp1*kp2*(2*(me*me) - 2*p1p3)*
             (2*(me*me) - 2*(kp1 - kp3 + p1p3)))) + 
       ((64*(2*pow(me,6) + pow(me,4)*
                (-kp1 - kp2 + kp3 - 7*(me*me) - 7*p1p2 - 
                  20*(-kp1 - kp2 + p1p2) + 
                  5*(-kp1 + me*me + p1p2 - p1p3) + 26*p1p3 + 
                  15*(kp1 - kp3 + p1p3) - 
                  2*(17*kp1 - 19*(kp1 + kp2 - kp3) - 11*kp3 - 
                     13*(me*me) - 18*(-kp1 - kp2 + p1p2) - 
                     7*(-kp1 + me*me + p1p2 - p1p3) + 17*p1p3)) - 
               4*(-kp1 + me*me + p1p2 - p1p3)*
                (-5*(me*me)*(-kp1 - kp2 + p1p2) - 
                  me*me*(kp2 + kp3 - 
                     6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) - 
                  10*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  2*(-kp1 - kp2 + p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  2*(-kp1 + me*me + p1p2 - p1p3)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                  me*me*(kp1 - kp3 + p1p3) + 
                  4*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                  2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (kp1 - kp3 + p1p3) + 
                  p1p2*(me*me - 4*(-kp1 - kp2 + p1p2) - 
                     2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     4*p1p3) + 
                  p1p3*(5*(me*me) - 
                     2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     4*(kp1 - kp3 + p1p3))) + 
               me*me*(-26*pow(me,4) + 6*(p1p2*p1p2) + 
                  14*(me*me)*(-kp1 - kp2 + p1p2) + 
                  4*pow(-kp1 - kp2 + p1p2,2) + 
                  27*(me*me)*(-kp1 + me*me + p1p2 - p1p3) + 
                  12*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) + 
                  31*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  12*(-kp1 - kp2 + p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  24*(-kp1 + me*me + p1p2 - p1p3)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  14*(p1p3*p1p3) + 29*(me*me)*(kp1 - kp3 + p1p3) + 
                  2*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                  38*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) - 
                  14*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (kp1 - kp3 + p1p3) - 14*pow(kp1 - kp3 + p1p3,2) + 
                  p1p3*(17*(me*me) + 10*(-kp1 - kp2 + p1p2) - 
                     30*(-kp1 + me*me + p1p2 - p1p3) - 
                     6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     20*(kp1 - kp3 + p1p3)) + 
                  2*p1p2*(9*(-kp1 - kp2 + p1p2) + 
                     9*(-kp1 + me*me + p1p2 - p1p3) - 
                     3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     4*(kp1 - kp3 + p1p3)))))/
           (kp2*kp2*(2*(me*me) - 2*p1p3)) - 
          (64*(4*pow(me,6) + 
               pow(me,4)*(4*(me*me) + 8*p1p2 + 
                  12*(-kp1 - kp2 + p1p2) - 
                  6*(-kp1 + me*me + p1p2 - p1p3) - 
                  3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 3*p1p3 - 
                  10*(kp1 - kp3 + p1p3) + 
                  2*(17*kp1 - 4*(kp1 + kp2 - kp3) - 8*(me*me) + 
                     5*(-kp1 + me*me + p1p2 - p1p3) + 19*p1p3)) + 
               me*me*(-10*pow(me,4) - 4*(p1p2*p1p2) + 
                  8*(me*me)*(-kp1 - kp2 + p1p2) + 
                  12*pow(-kp1 - kp2 + p1p2,2) + 
                  6*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
                  4*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) - 
                  33*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  30*(-kp1 - kp2 + p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  30*(-kp1 + me*me + p1p2 - p1p3)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  14*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                  2*(p1p3*p1p3) + 10*(me*me)*(kp1 - kp3 + p1p3) - 
                  4*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                  8*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                  14*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (kp1 - kp3 + p1p3) - 
                  4*(-2*kp1 - kp2 + kp3 + 2*(me*me) + 2*p1p2 - 2*p1p3)*
                   (-kp2 - kp3 + 2*p1p2 + 2*p1p3) - 
                  2*p1p2*(15*(me*me) - 4*(-kp1 - kp2 + p1p2) - 
                     8*(-kp1 + me*me + p1p2 - p1p3) - 
                     3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     11*p1p3 - 4*(kp1 - kp3 + p1p3)) - 
                  p1p3*(-21*(me*me) - 2*(-kp1 - kp2 + p1p2) + 
                     2*(-kp1 + me*me + p1p2 - p1p3) + 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     42*(kp1 - kp3 + p1p3))) + 
               4*(2*(-kp1 + me*me + p1p2 - p1p3)*
                   pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                  2*(kp2*kp2)*p1p3 + 
                  kp2*(me*me*(-kp1 + me*me) - (2*kp3 + me*me)*p1p2 - 
                     (2*kp3 + 3*(me*me))*p1p3) - 
                  (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (me*me*(-2*kp1 + kp2 + kp3 + 2*(me*me)) + 
                     2*(kp1 + kp3)*p1p2 - 2*(p1p2*p1p2) - 
                     2*(-kp1 + kp3 + 2*(me*me))*p1p3 + 2*(p1p3*p1p3)) + 
                  kp3*((2*kp3 + 3*(me*me))*p1p2 + 
                     me*me*(-kp1 + me*me + p1p3)))))/
           (kp1*kp2*(2*(me*me) - 2*(kp1 - kp3 + p1p3))) + 
          256*((28*pow(me,6) + 
                pow(me,4)*(6*(kp1 + kp2 - kp3) + kp3 + 3*(me*me) + 
                   3*(-kp1 - kp2 + p1p2) - 
                   2*(18*kp1 + 2*kp2 - 2*kp3 + 19*(me*me) + 3*p1p2 - 
                      32*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      6*p1p3) - 20*(-kp1 + me*me + p1p2 - p1p3) + 6*p1p3\
) - 8*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                 (kp1*kp1 + kp1*kp2 - kp1*kp3 - 2*kp1*(me*me) + 
                   (2*kp1 + kp2 - kp3)*(me*me) - pow(me,4) - 
                   p1p2*p1p2 - p1p2*(kp1 + me*me - p1p2 - p1p3) - 
                   kp1*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                   me*me*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   (-kp1 - me*me + p1p2 - p1p3)*p1p3 - p1p3*p1p3) + 
                me*me*(-14*pow(me,4) + 4*(p1p2*p1p2) + 
                   16*(me*me)*(-kp1 - kp2 + p1p2) + 
                   20*pow(-kp1 - kp2 + p1p2,2) + 
                   26*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
                   12*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3) - 
                   37*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                   34*(-kp1 - kp2 + p1p2)*
                    (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                   6*(-kp1 + me*me + p1p2 - p1p3)*
                    (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                   30*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                   2*(p1p3*p1p3) + 6*(me*me)*(kp1 - kp3 + p1p3) - 
                   28*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                   22*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (kp1 - kp3 + p1p3) + 8*pow(kp1 - kp3 + p1p3,2) + 
                   2*p1p2*(3*(me*me) + 8*(-kp1 - kp2 + p1p2) - 
                      2*(-kp1 + me*me + p1p2 - p1p3) - 
                      5*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      p1p3 - 2*(kp1 - kp3 + p1p3)) - 
                   p1p3*(7*(me*me) + 10*(-kp1 - kp2 + p1p2) + 
                      2*(-kp1 + me*me + p1p2 - p1p3) - 
                      4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      2*(kp1 - kp3 + p1p3))))/
              (4.*(kp1*kp1)*(2*(me*me) - 2*(kp1 - kp3 + p1p3))) + 
             (6*pow(me,6) - 
                pow(me,4)*(4*kp1 + 17*(kp1 + kp2 - kp3) - 5*kp3 + 
                   4*(me*me) + 4*(-kp1 + me*me + p1p2 - p1p3) + 4*p1p3 + 
                   2*(5*kp1 + 27*(kp1 + kp2 - kp3) - 5*kp3 + 7*(me*me) + 
                      4*(-kp1 - kp2 + p1p2) - 
                      9*(-kp1 + me*me + p1p2 - p1p3) + 11*p1p3)) + 
                me*me*(-4*(kp2*kp2) - 8*kp2*(kp1 + kp2 - kp3) - 
                   4*kp2*kp3 + 16*(kp1 + kp2 - kp3)*kp3 + 16*(kp3*kp3) - 
                   26*kp2*(me*me) + 67*(kp1 + kp2 - kp3)*(me*me) + 
                   25*kp3*(me*me) + 40*pow(me,4) - 
                   12*kp2*(-kp1 - kp2 + p1p2) - 
                   20*(kp1 + kp2 - kp3)*(-kp1 - kp2 + p1p2) + 
                   20*(me*me)*(-kp1 - kp2 + p1p2) - 
                   16*pow(-kp1 - kp2 + p1p2,2) + 
                   50*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                   (-38*kp2 + 22*(kp1 + kp2 - kp3) + 58*kp3 + 
                      24*(me*me) + 42*(-kp1 - kp2 + p1p2))*
                    (kp1 - kp3 + p1p3) - 26*pow(kp1 - kp3 + p1p3,2) - 
                   2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (-15*kp2 + 39*(kp1 + kp2 - kp3) + 29*kp3 + 
                      46*(me*me) + 41*(-kp1 - kp2 + p1p2) - 
                      4*(kp1 - kp3 + p1p3))) - 
                4*(2*(2*kp1 + me*me + p1p2 - p1p3)*
                    pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                   2*(kp2*kp2)*p1p3 + 
                   kp2*(me*me*(-kp1 + me*me) - (2*kp3 + me*me)*p1p2 - 
                      (2*kp3 + 3*(me*me))*p1p3) + 
                   2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                    (3*(kp1*kp1) + 3*kp1*kp2 - 3*kp1*kp3 - 
                      8*kp1*(me*me) - 2*kp2*(me*me) + kp3*(me*me) + 
                      2*pow(me,4) + (-4*kp1 - kp3 + 3*(me*me))*p1p2 + 
                      p1p2*p1p2 + (2*kp1 + kp3 - me*me)*p1p3 - p1p3*p1p3) \
+ me*me*(-((-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)) - 
                      4*(-kp1 - kp2 + p1p2)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      6*(-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      p1p2*(-kp1 + me*me + p1p2 - 
                         2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                         p1p3) - 
                      (-kp1 + me*me + p1p2 + 
                        4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        p1p3)*p1p3 - 
                      (-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                      2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                       (kp1 - kp3 + p1p3)) + 
                   kp3*((2*kp3 + 3*(me*me))*p1p2 + 
                      me*me*(-kp1 + me*me + p1p3))))/
              (4.*kp1*kp2*(2*(me*me) - 2*p1p3))))/
        (2.*(-kp1 - kp2 + me*me + p1p2)) + 
       ((-512*(-3*pow(me,6) - 2*pow(me,4)*(-kp3 + 2*p1p2) + 
               kp3*(p1p2*(-kp3 + me*me + p1p2) - 
                  2*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3)) - 
               me*me*(-2*kp3*p1p2 + p1p2*p1p2 - kp3*(kp1 - kp3 + p1p3) - 
                  (-kp1 + me*me + p1p2 - p1p3)*
                   (kp3 + 2*(kp1 - kp3 + p1p3)))))/pow(me*me + p1p2,2) + 
          (256*(12*pow(me,6) - 
               8*kp3*(-kp1 + kp3 + p1p2 - p1p3)*
                (-kp1 + me*me + p1p2 - p1p3) + 
               pow(me,4)*(30*kp2 + 4*(kp1 + kp2 - kp3) - 19*kp3 - 
                  13*(me*me) - 5*p1p2 + 8*(-kp1 - kp2 + p1p2) - 
                  3*(-kp1 + me*me + p1p2 - p1p3) + 
                  15*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  2*(kp1 - kp3 + p1p3)) + 
               me*me*(2*(kp1*kp1) + 4*kp1*kp2 + 2*(kp2*kp2) + 
                  4*kp1*(kp1 + kp2 - kp3) + 4*kp2*(kp1 + kp2 - kp3) - 
                  6*pow(kp1 + kp2 - kp3,2) + kp1*(me*me) - 
                  14*kp2*(me*me) - 5*(kp1 + kp2 - kp3)*(me*me) + 
                  14*pow(me,4) + 4*(p1p2*p1p2) + 
                  8*pow(-kp1 + me*me + p1p2 - p1p3,2) + 
                  10*kp1*(kp1 - kp3 + p1p3) + 
                  10*kp2*(kp1 - kp3 + p1p3) - 
                  2*(kp1 + kp2 - kp3)*(kp1 - kp3 + p1p3) - 
                  19*(me*me)*(kp1 - kp3 + p1p3) + 
                  4*pow(kp1 - kp3 + p1p3,2) + 
                  2*p1p2*(-4*kp1 - 4*kp2 - kp3 + 7*(me*me) - 
                     2*(-kp1 + me*me + p1p2 - p1p3) - 
                     4*(kp1 - kp3 + p1p3)) + 
                  (-kp1 + me*me + p1p2 - p1p3)*
                   (-2*kp1 - 2*kp2 + 10*(kp1 + kp2 - kp3) - 9*(me*me) + 
                     4*(kp1 - kp3 + p1p3)))))/
           ((me*me + p1p2)*(2*(me*me) - 2*(kp1 - kp3 + p1p3))) + 
          (128*(138*pow(me,6) + 
               pow(me,4)*(39*p1p2 + 74*(-kp1 - kp2 + p1p2) - 
                  63*(-kp1 + me*me + p1p2 - p1p3) - 
                  67*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                  37*p1p3 - 117*(kp1 - kp3 + p1p3)) + 
               me*me*(8*(kp1*kp1) + 8*kp1*kp2 + 
                  4*kp1*(kp1 + kp2 - kp3) + 12*kp2*(kp1 + kp2 - kp3) - 
                  12*pow(kp1 + kp2 - kp3,2) + 17*kp1*(me*me) - 
                  5*kp2*(me*me) + 10*pow(me,4) + 24*(p1p2*p1p2) + 
                  20*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
                  12*kp1*(kp1 - kp3 + p1p3) - 
                  20*kp2*(kp1 - kp3 + p1p3) + 
                  32*(kp1 + kp2 - kp3)*(kp1 - kp3 + p1p3) + 
                  26*(me*me)*(kp1 - kp3 + p1p3) - 
                  20*pow(kp1 - kp3 + p1p3,2) + 
                  2*(-kp1 + me*me + p1p2 - p1p3)*
                   (-2*kp1 - 6*kp2 + 12*(kp1 + kp2 - kp3) + 5*(me*me) - 
                     8*(kp1 - kp3 + p1p3)) + 
                  p1p2*(-16*kp1 - 8*kp2 - 4*(kp1 + kp2 - kp3) + 
                     11*(me*me) - 12*(-kp1 + me*me + p1p2 - p1p3) + 
                     12*(kp1 - kp3 + p1p3)) - 
                  4*(p1p2*(-7*(-kp1 - kp2 + p1p2) - 
                        4*(-kp1 + me*me + p1p2 - p1p3) + 
                        4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        4*p1p3 - 4*(kp1 - kp3 + p1p3)) - 
                     2*((-kp1 + me*me + p1p2 - p1p3)*
                         (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        p1p3*(kp1 - kp3 + p1p3)))) + 
               2*(p1p2*p1p2*(8*(me*me) - 7*(-kp1 - kp2 + p1p2) + 
                     8*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 8*p1p3) \
+ p1p2*(-44*(me*me)*(-kp1 - kp2 + p1p2) - 
                     15*pow(-kp1 - kp2 + p1p2,2) + 
                     23*(-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3) - 
                     16*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                       2) - 16*(p1p3*p1p3) - 
                     8*(me*me)*(kp1 - kp3 + p1p3) + 
                     15*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                     p1p3*(40*(me*me) + 31*(-kp1 - kp2 + p1p2) - 
                        24*(-kp1 + me*me + p1p2 - p1p3) - 
                        32*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        18*(kp1 - kp3 + p1p3)) + 
                     (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (40*(me*me) + 31*(-kp1 - kp2 + p1p2) - 
                        26*(-kp1 + me*me + p1p2 - p1p3) - 
                        16*(kp1 - kp3 + p1p3))) + 
                  2*(pow(-kp1 + me*me + p1p2 - p1p3,2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     p1p3*(kp1 - kp3 + p1p3)*(kp1 - 3*(me*me) + p1p3) + 
                     (-kp1 + me*me + p1p2 - p1p3)*
                      (pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                        (kp1 - kp3 + p1p3)*
                         (-4*(-kp1 - kp2 + me*me + p1p2) + 5*p1p3) + 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                         (kp1 + kp2 - 4*(me*me) - p1p2 + p1p3 + 
                         5*(kp1 - kp3 + p1p3)))))))/
           pow(2*(me*me) - 2*(kp1 - kp3 + p1p3),2))/(4.*(kp3*kp3)) + 
       ((128*((9*pow(me,6) + 
                  pow(me,4)*
                   (-44*kp1 - 21*kp2 - 39*(kp1 + kp2 - kp3) + 5*kp3 - 
                     54*(me*me) + 4*p1p2 - 53*(-kp1 - kp2 + p1p2) + 
                     16*(-kp1 + me*me + p1p2 - p1p3) + 
                     24*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     4*p1p3) + 
                  p1p2*p1p2*(8*kp3 - 10*(me*me) + 
                     26*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     10*p1p3) + 
                  p1p2*(16*kp2*kp3 - 8*(kp3*kp3) + 21*kp2*(me*me) + 
                     10*kp3*(me*me) + 
                     kp1*(-8*kp3 + 5*(me*me) - 
                        34*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) \
- 58*kp2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     34*kp3*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     26*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     26*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                       2) + 2*
                      (-9*kp1 + 13*(-kp2 + kp3 + me*me) - 
                       42*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                      p1p3 - 26*(p1p3*p1p3)) - 
                  me*me*(50*pow(me,4) + 14*(p1p2*p1p2) + 
                     46*(me*me)*(-kp1 - kp2 + p1p2) + 
                     pow(-kp1 - kp2 + p1p2,2) - 
                     42*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
                     (-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3) - 
                     77*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     29*(-kp1 - kp2 + p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     52*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     12*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                       2) + p1p3*p1p3 - 
                     43*(me*me)*(kp1 - kp3 + p1p3) - 
                     6*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                     5*(-kp1 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) + 
                     25*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) - 
                     3*pow(kp1 - kp3 + p1p3,2) - 
                     8*((-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                       (-kp1 + kp3 + 3*(-kp1 - kp2 + p1p2) + 
                       2*(-kp1 + me*me + p1p2 - p1p3) - p1p3) + 
                        (kp2 + me*me + 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                       p1p3)*p1p3) + 
                     p1p3*(-20*(me*me) - 10*(-kp1 - kp2 + p1p2) + 
                        9*(-kp1 + me*me + p1p2 - p1p3) + 
                        37*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        34*(kp1 - kp3 + p1p3)) + 
                     p1p2*(70*(me*me) + 55*(-kp1 - kp2 + p1p2) - 
                        30*(-kp1 + me*me + p1p2 - p1p3) - 
                        42*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        39*p1p3 - 27*(kp1 - kp3 + p1p3)) + 
                     p1p2*(-10*(-kp1 - kp2 + p1p2) - 
                        21*(-kp1 + me*me + p1p2 - p1p3) + 
                        13*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        13*p1p3 - 5*(kp1 - kp3 + p1p3))) + 
                  8*(((2*kp2 - kp3)*(-kp2 + me*me) + 
                        kp1*
                       (-kp2 + me*me + 
                       4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                        (4*kp2 - 3*kp3 - 4*(me*me))*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2))*p1p3 + 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (p1p3*p1p3) + 
                     (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (kp1*kp1 + 
                        me*me*
                        (-3*kp2 + kp3 - 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                        kp1*(kp1 + 3*kp2 + me*me - p1p2 + p1p3))))/
                pow(2*(me*me) - 2*p1p3,2) + 
               (6*pow(me,6) + 
                  pow(me,4)*
                   (-3*kp1 + 14*kp2 - 21*(me*me) - 17*p1p2 + 
                     12*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     10*p1p3 - 
                     2*(kp1 + kp2 + 13*(kp1 + kp2 - kp3) - 7*kp3 + 
                        2*(me*me) - p1p2 + 
                        5*(-kp1 + me*me + p1p2 - p1p3) + 3*p1p3)) + 
                  8*(p1p2*((kp2 - kp3)*(kp3 + me*me) + 
                        (-2*kp1 - 4*kp2 + 3*kp3 + 2*(me*me) - 
                       6*p1p3)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        p1p3,2)) + 
                     2*(p1p2*p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     ((kp2 - kp3)*(-kp2 + me*me) + 
                        (4*kp1 + 4*kp2 - 3*kp3 - 4*(me*me))*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2))*p1p3 + 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (p1p3*p1p3) + 
                     (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (me*me*
                        (-2*kp2 + kp3 - 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                        kp1*(kp1 + 2*kp2 + 2*(me*me) - p1p2 + p1p3))) + 
                  me*me*(-4*pow(me,4) - 2*(p1p2*p1p2) + 
                     46*(me*me)*(-kp1 - kp2 + p1p2) + 
                     48*pow(-kp1 - kp2 + p1p2,2) + 
                     45*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
                     24*(-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3) - 
                     24*pow(-kp1 + me*me + p1p2 - p1p3,2) + 
                     4*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     18*(-kp1 - kp2 + p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     14*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     10*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                       2) + 8*(p1p3*p1p3) + 
                     12*(me*me)*(kp1 - kp3 + p1p3) - 
                     46*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                     18*(-kp1 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) + 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) + 
                     6*pow(kp1 - kp3 + p1p3,2) - 
                     2*p1p2*(kp1 - p1p2 - 3*(-kp1 - kp2 + p1p2) + 
                        4*p1p3) - 
                     p1p3*(17*(me*me) + 44*(-kp1 - kp2 + p1p2) + 
                        12*(-kp1 + me*me + p1p2 - p1p3) + 
                        6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        30*(kp1 - kp3 + p1p3)) - 
                     8*((-kp2 + kp3)*p1p2 - 
                        (-kp1 + me*me + p1p2 + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        p1p3)*p1p3 + 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                         (kp1 - kp3 - 2*(-kp1 - kp2 + p1p2) - 
                         2*(-kp1 + me*me + p1p2 - p1p3) + p1p3))))/
                (2.*(me*me + p1p2)*(2*(me*me) - 2*p1p3))))/kp2 - 
          ((-128*(20*pow(me,6) - 
                  4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                   (-2*(p1p2*p1p2) + 
                     me*me*(2*kp1 + kp2 - kp3 + 
                        2*(-kp1 + me*me + p1p2 - p1p3)) - 
                     me*me*(2*kp1 - kp2 + kp3 + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) - 
                     2*p1p2*
                      (-2*kp1 - kp2 + kp3 + 2*(me*me) + p1p2 - 
                       3*p1p3) + 
                     2*(-2*kp1 - kp2 + kp3 + 2*(me*me) + p1p2 - 
                       p1p3)*p1p3 - 2*(p1p3*p1p3)) + 
                  pow(me,4)*
                   (8*(me*me) + 8*p1p2 - 
                     13*(-kp1 + me*me + p1p2 - p1p3) - 
                     9*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     p1p3 - 13*(kp1 - kp3 + p1p3) + 
                     2*(2*(me*me) + 12*p1p2 + 6*(-kp1 - kp2 + p1p2) - 
                        5*(-kp1 + me*me + p1p2 - p1p3) + 
                        12*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        5*(kp1 - kp3 + p1p3))) + 
                  me*me*(4*pow(me,4) + 4*(p1p2*p1p2) + 
                     18*(me*me)*(-kp1 - kp2 + p1p2) + 
                     22*pow(-kp1 - kp2 + p1p2,2) + 
                     9*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
                     10*(-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3) - 
                     4*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
                     57*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     36*(-kp1 - kp2 + p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     46*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     26*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                       2) + 14*(p1p3*p1p3) + 
                     13*(me*me)*(kp1 - kp3 + p1p3) - 
                     6*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                     12*(-kp1 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) + 
                     22*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) - 
                     8*pow(kp1 - kp3 + p1p3,2) - 
                     p1p3*(21*(me*me) + 12*(-kp1 - kp2 + p1p2) + 
                        2*(-kp1 + me*me + p1p2 - p1p3) - 
                        8*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        10*(kp1 - kp3 + p1p3)) - 
                     2*p1p2*(5*p1p3 + 
                        3*(-2*kp1 - 2*kp2 + kp3 + 2*(me*me) + 2*p1p2 - 
                        2*(-kp1 + me*me + p1p2 - p1p3) - p1p3 - 
                        2*(kp1 - kp3 + p1p3))))))/(me*me + p1p2) - 
             (128*(30*pow(me,6) + 
                  pow(me,4)*
                   (45*(me*me) + 52*p1p2 + 25*(-kp1 - kp2 + p1p2) - 
                     5*(-kp1 + me*me + p1p2 - p1p3) + 
                     26*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     60*p1p3 - 47*(kp1 - kp3 + p1p3)) + 
                  4*(2*(2*(kp1*kp1) + 2*kp1*kp2 + 
                        me*me*
                        (-3*kp2 + kp3 + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)))*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     2*(p1p2*p1p2)*
                      (2*kp3 + me*me - 
                        3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        p1p3) + 
                     4*(kp2*(-kp2 + me*me) - 
                        me*me*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2) + 
                        kp1*
                        (-kp1 - 2*kp2 + kp3 + 2*(me*me) + p1p2 - 
                       p1p3))*p1p3 + 
                     4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (p1p3*p1p3) + 
                     p1p2*(4*kp2*kp3 + kp2*(me*me) - 2*kp3*(me*me) + 
                        kp1*(-4*kp3 - 3*(me*me) + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                        2*kp2*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        6*kp3*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        10*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        p1p3,2) + 
                        2*(3*kp1 + kp2 + kp3 + me*me - 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        p1p3 - 2*(p1p3*p1p3))) + 
                  me*me*(4*(kp2*kp2) - 4*pow(kp1 + kp2 - kp3,2) + 
                     4*kp2*kp3 - 28*(kp1 + kp2 - kp3)*kp3 - 
                     8*(kp3*kp3) - 4*kp2*(me*me) - 
                     19*(kp1 + kp2 - kp3)*(me*me) - 10*pow(me,4) + 
                     4*kp2*(-kp1 - kp2 + p1p2) - 
                     60*(kp1 + kp2 - kp3)*(-kp1 - kp2 + p1p2) - 
                     4*(me*me)*(-kp1 - kp2 + p1p2) + 
                     8*pow(-kp1 - kp2 + p1p2,2) + 
                     54*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                       2) + 2*
                      (kp2 - 13*(kp1 + kp2 - kp3) + 5*kp3 + 
                        22*(me*me) + 13*(-kp1 - kp2 + p1p2))*
                      (kp1 - kp3 + p1p3) - 
                     2*pow(kp1 - kp3 + p1p3,2) - 
                     2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (-5*kp2 - 23*(kp1 + kp2 - kp3) + 15*kp3 + 
                        22*(me*me) + 31*(-kp1 - kp2 + p1p2) + 
                        22*(kp1 - kp3 + p1p3)) - 
                     4*(-2*(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        2*(-kp1 - kp2 + p1p2) - p1p3)*p1p3 + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 - 3*(-kp1 - kp2 + p1p2) + 
                        2*(-kp1 + me*me + p1p2 - p1p3) + p1p3) + 
                        p1p2*
                         (2*kp1 + kp2 - kp3 - 2*(me*me) - 2*p1p2 - 
                         2*(-kp1 - kp2 + p1p2) + 3*p1p3 + 
                         3*(kp1 - kp3 + p1p3))))))/(2*(me*me) - 2*p1p3))/
           (2.*kp1*(2*(me*me) - 2*(kp1 - kp3 + p1p3))) + 
          ((-128*(-18*pow(me,6) + 
                  pow(me,4)*
                   (29*(me*me) + 68*p1p2 + 53*(-kp1 - kp2 + p1p2) - 
                     37*(-kp1 + me*me + p1p2 - p1p3) - 
                     10*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     16*p1p3 + 5*(kp1 - kp3 + p1p3)) + 
                  16*(kp3*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     2*(p1p2*p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     kp1*(-kp1 + 2*kp3 + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     kp2*kp2*p1p3 - kp2*kp3*p1p3 - kp2*(me*me)*p1p3 + 
                     kp3*(me*me)*p1p3 + 
                     kp2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      p1p3 - 
                     2*kp3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      p1p3 + 
                     p1p2*(-(kp2*kp3) + kp3*kp3 - 2*kp2*(me*me) - 
                        kp3*(me*me) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2) - 
                        kp1*
                        (-kp1 - kp2 + kp3 + 2*(me*me) + p1p2 - p1p3) \
+ (kp1 + kp2 - 2*kp3 - 2*(me*me) + 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        p1p3 + 2*(p1p3*p1p3))) + 
                  me*me*(-4*kp1*kp2 - 20*(kp2*kp2) + 4*kp1*kp3 + 
                     24*kp2*kp3 - 4*(kp3*kp3) + 29*kp1*(me*me) + 
                     27*kp2*(me*me) - kp3*(me*me) - 8*pow(me,4) - 
                     20*(p1p2*p1p2) - 
                     34*kp1*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     6*kp2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     22*kp3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     30*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     14*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                       2) - (30*kp1 - 46*kp2 + 14*kp3 - 22*(me*me))*
                      p1p3 - 46*(p1p3*p1p3) - 
                     p1p2*(-20*kp1 - 40*kp2 + 24*kp3 + 42*(me*me) - 
                        6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        30*p1p3) + 
                     16*(p1p2*
                        (-kp1 + kp3 - 2*(-kp1 + me*me + p1p2 - p1p3)) \
- (-kp1 + me*me + p1p2 - p1p3)*p1p3 - 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3)))))/
              ((2*(me*me) - 2*p1p3)*(2*(me*me) - 2*(kp1 - kp3 + p1p3))) \
- (512*(2*pow(me,6) + 2*(p1p2*p1p2)*
                   (-kp1 - kp2 + kp3 + me*me + p1p2) + 
                  pow(me,4)*
                   (2*kp1 + 2*kp2 + 2*kp3 + 8*(-kp1 - kp2 + p1p2)) + 
                  (kp1*kp1 + 2*kp3*(me*me) - 
                     kp1*(-kp1 + 3*kp3 + me*me + p1p2 - p1p3))*
                   (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                  (kp2*kp2 - 2*kp2*kp3 + 2*kp3*(me*me) + 
                     kp1*(-kp1 - 2*kp2 + kp3 + me*me + p1p2 - 
                      p1p3) + 
                     (kp2 - 4*kp3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*p1p3 - 
                  kp2*(p1p3*p1p3) + 
                  p1p2*(-(kp2*kp3) + 2*(kp3*kp3) - kp2*(me*me) + 
                     2*kp3*(me*me) - 
                     kp1*(kp3 + me*me + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) - 
                     kp2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     2*kp3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     2*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     (kp1 - 
                        2*(-kp2 + kp3 + me*me - 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)))*
                      p1p3) + 
                  me*me*(-(kp1*kp3) - kp2*kp3 + 2*(kp3*kp3) - 
                     kp1*(me*me) - kp2*(me*me) + 2*kp3*(me*me) + 
                     2*(p1p2*p1p2) - 
                     p1p2*(-kp1 - kp2 + 2*(me*me) + 2*p1p2) - 
                     p1p2*(kp1 + kp2 + 4*kp3 + 4*(me*me) - 
                        4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        4*p1p3) + 
                     kp1*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     2*kp2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     2*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     (-2*kp1 + kp2 - 2*(me*me) + 
                        4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                      p1p3 - 
                     2*((-kp1 + me*me + p1p2 - p1p3)*p1p3 + 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3)))))/pow(me*me + p1p2,2) + 
             ((-256*(14*pow(me,6) + 
                     8*(p1p2*
                       (-((kp2 - kp3)*(kp3 + me*me)) + 
                       (-2*kp1 - kp2 + 2*kp3 + 2*(me*me) - 2*p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                        kp3*(me*me)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*(p1p2*p1p2)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        kp1*(-kp1 + 2*kp3 + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        kp2*kp2*p1p3 - kp2*kp3*p1p3 - 
                        kp2*(me*me)*p1p3 + kp3*(me*me)*p1p3 + 
                        kp2*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*p1p3 \
- 2*kp3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*p1p3) + 
                     pow(me,4)*
                      (kp2 - 12*kp3 - 7*(me*me) + 17*p1p2 + 
                        16*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*p1p3 + 
                        2*(4*kp1 - 2*(kp1 + kp2 - kp3) - 2*kp3 + 
                        2*(me*me) + 9*(-kp1 - kp2 + p1p2) - 
                        17*(-kp1 + me*me + p1p2 - p1p3) + 5*p1p3)) + 
                     me*me*(4*(kp1*kp1) + 4*(kp2*kp2) - 
                        8*kp1*(kp1 + kp2 - kp3) + 
                        8*pow(kp1 + kp2 - kp3,2) + 6*kp1*(me*me) - 
                        3*kp2*(me*me) - 4*(kp1 + kp2 - kp3)*(me*me) - 
                        3*pow(me,4) + 8*(p1p2*p1p2) - 
                        8*pow(-kp1 + me*me + p1p2 - p1p3,2) + 
                        28*(kp1 + kp2 - kp3)*(kp1 - kp3 + p1p3) - 
                        20*(me*me)*(kp1 - kp3 + p1p3) + 
                        4*pow(kp1 - kp3 + p1p3,2) + 
                        2*(-kp1 + me*me + p1p2 - p1p3)*
                        (-2*kp1 + 2*kp2 - 12*(kp1 + kp2 - kp3) + 
                        9*(me*me) + 6*(kp1 - kp3 + p1p3)) - 
                        p1p2*
                        (-12*kp1 + 12*kp2 + 32*(kp1 + kp2 - kp3) + 
                        19*(me*me) - 
                        48*(-kp1 + me*me + p1p2 - p1p3) + 
                        12*(kp1 - kp3 + p1p3)) - 
                        8*(p1p2*
                        (-2*kp1 - kp2 + kp3 + 2*(me*me) + 2*p1p2 - 
                        2*p1p3) + (-kp1 + me*me + p1p2 - p1p3)*p1p3 + 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3)))))/(2*(me*me) - 2*p1p3) - 
                (512*(4*pow(me,6) + 
                     pow(me,4)*
                      (-kp1 - kp2 + 17*kp3 + me*me + 
                        17*(-kp1 - kp2 + p1p2) + 
                        3*(-kp1 + me*me + p1p2 - p1p3) - 
                        16*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        p1p3) + 
                     4*(p1p2*
                        (-((kp2 - kp3)*(kp3 + me*me)) + 
                        (-2*kp1 - kp2 + 2*kp3 + 2*(me*me) - 2*p1p3)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                        kp3*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*(p1p2*p1p2)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        kp1*(-kp1 + 2*kp3 + p1p2 - p1p3)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        kp2*kp2*p1p3 - kp2*kp3*p1p3 - 
                        kp2*(me*me)*p1p3 + kp3*(me*me)*p1p3 + 
                        kp2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        p1p3 - 
                        2*kp3*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*p1p3) \
+ me*me*(4*(kp1*kp1) + 8*kp2*(kp1 + kp2 - kp3) - 
                        2*pow(kp1 + kp2 - kp3,2) + 2*kp1*(me*me) + 
                        3*kp2*(me*me) - 15*(kp1 + kp2 - kp3)*(me*me) - 
                        2*pow(me,4) + 4*(p1p2*p1p2) - 
                        2*p1p2*
                        (11*(kp1 + kp2 - kp3) + 5*(me*me) - 
                        14*(-kp1 + me*me + p1p2 - p1p3)) - 
                        8*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
                        4*kp2*(kp1 - kp3 + p1p3) + 
                        18*(kp1 + kp2 - kp3)*(kp1 - kp3 + p1p3) + 
                        2*(me*me)*(kp1 - kp3 + p1p3) - 
                        4*pow(kp1 - kp3 + p1p3,2) + 
                        2*(-kp1 + me*me + p1p2 - p1p3)*
                        (-4*kp1 - 4*kp2 + 3*(kp1 + kp2 - kp3) + 
                        9*(me*me) - 2*(kp1 - kp3 + p1p3)) - 
                        4*(p1p2*
                        (-2*kp1 - kp2 + kp3 + 2*(me*me) + 2*p1p2 - 
                        2*p1p3) + (-kp1 + me*me + p1p2 - p1p3)*p1p3 + 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3)))))/
                 (2*(me*me) - 2*(kp1 - kp3 + p1p3)))/(2.*(me*me + p1p2)))/
           (2.*kp3) + (-((-512*
                    (13*pow(me,6) - 
                      2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                       (-2*(p1p2*p1p2) + 
                        me*me*
                       (2*kp1 + kp2 - kp3 + 
                       2*(-kp1 + me*me + p1p2 - p1p3)) - 
                        me*me*
                       (2*kp1 - kp2 + kp3 + 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) - 
                        2*p1p2*
                       (-2*kp1 - kp2 + kp3 + 2*(me*me) + p1p2 - 
                       3*p1p3) + 
                        2*(-2*kp1 - kp2 + kp3 + 2*(me*me) + p1p2 - 
                       p1p3)*p1p3 - 2*(p1p3*p1p3)) + 
                      me*me*
                       (4*(kp1*kp1) - kp1*kp2 - kp2*kp2 - 
                        4*kp1*kp3 + 5*kp2*kp3 + 3*kp2*(me*me) - 
                        2*kp3*(me*me) - pow(me,4) + 
                        8*(p1p2*p1p2) + 
                        p1p2*
                       (-kp2 - 2*kp3 + me*me + 
                       14*(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                      p1p3) - 14*p1p3) - 
                        10*kp1*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        3*kp2*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        4*kp3*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        14*(me*me)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2) + 
                        (6*kp1 + 3*kp2 - 10*(me*me) - 
                       16*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                       p1p3 + 6*(p1p3*p1p3)) + 
                      pow(me,4)*
                       (9*kp1 + 3*kp2 - 11*(kp1 + kp2 - kp3) - 
                        17*(me*me) + 
                        11*(-kp1 + me*me + p1p2 - p1p3) + 
                        15*(kp1 - kp3 + p1p3) + 
                        2*(kp1 - 4*kp2 - me*me + 
                        2*(-kp1 + me*me + p1p2 - p1p3) - 
                        2*(kp1 - kp3 + p1p3)))))/(2*(me*me) - 2*p1p3) - 
                 (512*(8*pow(me,6) + 
                      2*(p1p2*p1p2)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2) + 
                      pow(me,4)*
                       (-kp1 - kp2 + 2*kp3 + 2*p1p2 + 
                        2*(-kp1 - kp2 + p1p2) + 
                        2*(kp1 + 2*(me*me) + 5*p1p2 + 
                        3*(-kp1 - kp2 + p1p2) - 2*p1p3) - 4*p1p3) + 
                      (2*(kp1*kp1) - 2*kp1*kp3 + 
                        me*me*
                        (-kp2 + 3*kp3 + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)))*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      (kp1*(-2*kp2 + me*me) + kp3*(-2*kp2 + me*me) - 
                        2*kp2*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        4*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2))*p1p3 + 
                      (2*kp2 - 2*(me*me) + 
                        4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                       (p1p3*p1p3) + 
                      2*p1p2*
                       (kp3*kp3 + kp3*(me*me) + 
                        me*me*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2) - 
                        kp1*
                        (-kp1 - kp2 + 2*kp3 + 2*(me*me) + p1p2 - 
                       p1p3) + 
                        (kp1 - kp2 + me*me - 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        p1p3 - p1p3*p1p3) + 
                      me*me*(-2*kp1*kp3 + 2*(kp3*kp3) - kp1*(me*me) - 
                        kp2*(me*me) + 2*kp3*(me*me) + 4*(p1p2*p1p2) + 
                        2*kp1*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*kp2*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        4*kp3*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        p1p3,2) - 
                        2*p1p2*
                        (-2*kp1 + kp2 + 4*kp3 + 2*(me*me) - 
                        3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        p1p3) + (2*kp1 + 4*kp2 - 2*(me*me))*p1p3 - 
                        2*(p1p3*p1p3) - 
                        p1p3*
                        (-2*kp1 - kp2 + me*me + 2*p1p2 - p1p3 - 
                        2*(kp1 - kp3 + p1p3)) - 
                        p1p2*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 + 
                        2*(kp1 - kp3 + p1p3)) - 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 + kp2 - p1p2 + 
                        2*(-kp1 + me*me + p1p2 - p1p3) + 
                        3*(kp1 - kp3 + p1p3)))))/(me*me + p1p2))/(2.*kp1) \
- ((512*(8*pow(me,6) + 2*(p1p2*p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2) + 
                     ((kp2 + kp3)*(-2*kp1 + me*me) - 
                        2*(-kp1 + me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     (2*(kp2*kp2) - 2*kp2*kp3 + 3*kp3*(me*me) - 
                        kp1*
                       (me*me + 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                        4*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2))*p1p3 + 
                     2*(me*me - 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                      (p1p3*p1p3) + 
                     pow(me,4)*
                      (kp3 - me*me + p1p2 + 2*(-kp1 - kp2 + p1p2) - 
                        3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*(-kp1 + kp3 + 3*(me*me) + 6*p1p2 + 
                        3*(-kp1 - kp2 + p1p2) - 
                        3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        p1p3) + p1p3) + 
                     2*p1p2*(-((kp2 - kp3)*(kp3 + me*me)) + 
                        (-kp1 + kp2 + me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2) + 
                        (-kp2 + me*me - 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        p1p3 + p1p3*p1p3) + 
                     me*me*(-2*kp2*kp3 + 2*(kp3*kp3) - kp1*(me*me) - 
                        kp2*(me*me) + 2*kp3*(me*me) + 4*(p1p2*p1p2) - 
                        p1p2*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 + 
                        2*(-kp1 + me*me + p1p2 - p1p3)) - 
                        2*p1p2*
                        (2*kp1 - kp2 + 3*kp3 + me*me - p1p2 - 2*p1p3) \
+ 4*kp1*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*kp2*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        p1p3,2) - 
                        2*(-kp1 - kp2 + 2*kp3 + me*me)*p1p3 + 
                        2*(p1p3*p1p3) - 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (-kp2 - kp3 + p1p2 - 
                        2*(-kp1 + me*me + p1p2 - p1p3) + p1p3) - 
                        p1p3*
                        (kp1 + kp2 - p1p2 + 
                        3*(-kp1 + me*me + p1p2 - p1p3) + 
                        2*(kp1 - kp3 + p1p3)))))/(me*me + p1p2) + 
                (512*(9*pow(me,6) + 
                     2*(2*p1p2*
                        (-((kp2 - kp3)*(kp3 + me*me)) + 
                        (-kp1 + kp2 + me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2)) + 
                        2*(p1p2*p1p2)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        ((kp2 + kp3)*(-2*kp1 + me*me) - 
                        2*(-kp1 + me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*(-((kp2 - kp3)*(-kp2 + me*me)) + 
                        (-kp1 - kp2 + me*me)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2))*p1p3 - 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (p1p3*p1p3)) + 
                     pow(me,4)*
                      (-kp2 - kp3 - 4*(me*me) + 15*p1p2 + 
                        5*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        3*p1p3 + 2*(kp1 - kp3 + p1p3)) + 
                     me*me*(-(kp1*kp2) - 5*(kp2*kp2) + 5*kp2*kp3 + 
                        3*kp1*(me*me) + kp2*(me*me) - 2*kp3*(me*me) + 
                        12*kp1*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        5*kp2*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        4*kp3*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        8*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2) + 
                        (-8*kp1 + 5*kp2 + 2*(me*me) + 
                        16*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        p1p3 - 6*(p1p3*p1p3) + 
                        p1p2*
                        (-2*kp1 + 9*kp2 - 18*kp3 - 10*(me*me) + 
                        10*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        6*p1p3) - 
                        2*(p1p2*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 + 
                        2*(-kp1 + me*me + p1p2 - p1p3) - p1p3) + 
                         (-kp1 - kp2 + kp3 + me*me + p1p2 + 
                        2*(-kp1 + me*me + p1p2 - p1p3) - p1p3)*p1p3 + 
                         (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                         (-kp2 - kp3 + p1p2 - 
                         2*(-kp1 + me*me + p1p2 - p1p3) + p1p3)))))/
                 (2*(me*me) - 2*p1p3))/(2.*kp2))/
           (2.*(-kp1 - kp2 + me*me + p1p2)))/(2.*(kp1 + kp2 - kp3)) + 
       (-((128*(-2*pow(me,6) + 
                   pow(me,4)*
                    (-8*kp1 + 24*kp3 + 18*(me*me) + 
                      38*(-kp1 - kp2 + p1p2) - 
                      15*(-kp1 + me*me + p1p2 - p1p3) - 29*p1p3 + 
                      2*(12*kp1 - 33*(kp1 + kp2 - kp3) - 13*kp3 - 
                        18*(me*me) - 23*(-kp1 - kp2 + p1p2) + 
                        7*(-kp1 + me*me + p1p2 - p1p3) + 9*p1p3)) - 
                   8*(-kp1 + me*me + p1p2 - p1p3)*
                    (me*me*(kp2 + kp3 - 
                        4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                      2*(-kp1 - kp2 + p1p2)*
                       (-kp1 + kp3 + me*me + p1p2 - p1p3) - 
                      2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2) - 
                      p1p3*(2*(me*me) + 2*p1p2 + 
                        3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*(kp1 - kp3 + p1p3)) - 
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                       (-10*(me*me) - 3*p1p2 - 3*(-kp1 - kp2 + p1p2) + 
                        4*(-kp1 + me*me + p1p2 - p1p3) + 
                        3*(kp1 - kp3 + p1p3))) + 
                   me*me*(-10*(kp1*kp1) - 24*kp1*kp2 - 14*(kp2*kp2) + 
                      6*kp1*(kp1 + kp2 - kp3) + 
                      6*kp2*(kp1 + kp2 - kp3) - 
                      4*pow(kp1 + kp2 - kp3,2) + 5*kp1*(me*me) + 
                      52*kp2*(me*me) + 12*(kp1 + kp2 - kp3)*(me*me) - 
                      5*pow(me,4) - 
                      60*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
                      10*kp1*(kp1 - kp3 + p1p3) - 
                      22*kp2*(kp1 - kp3 + p1p3) + 
                      44*(me*me)*(kp1 - kp3 + p1p3) - 
                      8*pow(kp1 - kp3 + p1p3,2) + 
                      (-kp1 + me*me + p1p2 - p1p3)*
                       (-54*kp1 - 34*kp2 + 52*(kp1 + kp2 - kp3) + 
                        38*(me*me) - 20*(kp1 - kp3 + p1p3)) + 
                      p1p2*(10*kp1 + 30*kp2 - 17*(me*me) + 
                         12*(-kp1 + me*me + p1p2 - p1p3) + 
                         8*(kp1 - kp3 + p1p3)))))/(me*me + p1p2) - 
              (128*(110*pow(me,6) + 
                   pow(me,4)*
                    (16*kp1 + 76*kp2 - 54*kp3 - 91*(me*me) + 28*p1p2 + 
                      13*(-kp1 - kp2 + p1p2) + 
                      10*(-kp1 + me*me + p1p2 - p1p3) + 
                      131*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      85*p1p3 - 82*(kp1 - kp3 + p1p3)) + 
                   2*(10*(p1p2*p1p2)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        3*(-kp1 - kp2 + p1p2)) + 
                      24*(me*me)*(-kp1 - kp2 + p1p2)*
                       (-kp1 + me*me + p1p2 - p1p3) - 
                      8*(me*me)*(-kp1 - kp2 + p1p2)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      76*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      27*(-kp1 - kp2 + p1p2)*
                       (-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      27*pow(-kp1 + me*me + p1p2 - p1p3,2)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      11*(-kp1 + me*me + p1p2 - p1p3)*
                       pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) \
- 8*(me*me)*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                      16*(-kp1 - kp2 + p1p2)*
                       (-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) \
- 19*(-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                       (kp1 - kp3 + p1p3) - 
                      13*(p1p3*p1p3)*(kp1 - kp3 + p1p3) + 
                      p1p2*(-118*(me*me)*(-kp1 - kp2 + p1p2) - 
                        38*pow(-kp1 - kp2 + p1p2,2) + 
                        (-13*(me*me) + 54*(-kp1 - kp2 + p1p2) + 
                        9*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        (-kp1 + me*me + p1p2 - p1p3) + 
                        25*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        48*(-kp1 - kp2 + p1p2)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        10*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        p1p3,2) - 18*(p1p3*p1p3) - 
                        5*(me*me)*(kp1 - kp3 + p1p3) + 
                        38*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                        10*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) + 
                        p1p3*
                        (49*(me*me) + 64*(-kp1 - kp2 + p1p2) - 
                        34*(-kp1 + me*me + p1p2 - p1p3) - 
                        28*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        5*(kp1 - kp3 + p1p3))) + 
                      p1p3*(8*(me*me)*(-kp1 - kp2 + p1p2) + 
                         (36*(me*me) + 5*(-kp1 - kp2 + p1p2) - 
                        13*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        (kp1 - kp3 + p1p3) - 
                         13*pow(kp1 - kp3 + p1p3,2) + 
                         (-kp1 + me*me + p1p2 - p1p3)*
                         (-16*(me*me) - 
                         19*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                         3*(kp1 - kp3 + p1p3)))) + 
                   me*me*(-117*pow(me,4) + 98*(p1p2*p1p2) - 
                      67*(me*me)*(-kp1 - kp2 + p1p2) + 
                      82*(me*me)*(-kp1 + me*me + p1p2 - p1p3) + 
                      38*(-kp1 - kp2 + p1p2)*
                       (-kp1 + me*me + p1p2 - p1p3) - 
                      6*pow(-kp1 + me*me + p1p2 - p1p3,2) + 
                      119*(me*me)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      44*(-kp1 - kp2 + p1p2)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      98*(-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      44*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2) - 28*(p1p3*p1p3) + 
                      66*(me*me)*(kp1 - kp3 + p1p3) - 
                      6*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                      16*(-kp1 + me*me + p1p2 - p1p3)*
                       (kp1 - kp3 + p1p3) - 
                      38*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                       (kp1 - kp3 + p1p3) + 
                      22*pow(kp1 - kp3 + p1p3,2) + 
                      2*p1p2*
                       (30*(me*me) + 17*(-kp1 - kp2 + p1p2) - 
                         14*(-kp1 + me*me + p1p2 - p1p3) - 
                         19*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                         35*p1p3 - 52*(kp1 - kp3 + p1p3)) + 
                      p1p3*(163*(me*me) + 44*(-kp1 - kp2 + p1p2) - 
                         66*(-kp1 + me*me + p1p2 - p1p3) - 
                         88*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                         42*(kp1 - kp3 + p1p3)) - 
                      2*(8*(-kp1 - kp2 + p1p2)*
                         (-kp1 + me*me + p1p2 - p1p3) - 
                         8*(-kp1 - kp2 + p1p2)*
                         (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                         38*(-kp1 + me*me + p1p2 - p1p3)*
                         (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                         8*(-kp1 + me*me + p1p2 - p1p3)*
                         (kp1 - kp3 + p1p3) + 
                         p1p2*
                         (-50*(-kp1 - kp2 + p1p2) - 
                         21*(-kp1 + me*me + p1p2 - p1p3) + 
                         5*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                         13*p1p3 - 5*(kp1 - kp3 + p1p3)) + 
                         p1p3*
                         (8*(-kp1 - kp2 + p1p2) - 
                         8*(-kp1 + me*me + p1p2 - p1p3) + 
                         10*(kp1 - kp3 + p1p3))))))/
               (2*(me*me) - 2*(kp1 - kp3 + p1p3)))/
           (2.*kp2*(2*(me*me) - 2*p1p3)) + 
          (-((-256*(12*pow(me,6) - 
                      4*(-kp1 + me*me + p1p2 - p1p3)*
                       (-(me*me*(-kp1 - kp2 + p1p2)) + 
                        me*me*
                       (kp2 + kp3 + 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) - 
                        2*(me*me)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*(-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2) + 3*(me*me)*(kp1 - kp3 + p1p3) + 
                        4*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                        p1p2*
                       (-3*(me*me) - 4*(-kp1 - kp2 + p1p2) + 4*p1p3) \
+ p1p3*(me*me - 4*(kp1 - kp3 + p1p3))) + 
                      me*me*
                       (22*pow(me,4) + 12*(p1p2*p1p2) + 
                        8*(me*me)*(-kp1 - kp2 + p1p2) + 
                        4*pow(-kp1 - kp2 + p1p2,2) - 
                        61*(me*me)*(-kp1 + me*me + p1p2 - p1p3) - 
                        18*(-kp1 - kp2 + p1p2)*
                       (-kp1 + me*me + p1p2 - p1p3) + 
                        38*pow(-kp1 + me*me + p1p2 - p1p3,2) + 
                        10*(me*me)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        4*(-kp1 - kp2 + p1p2)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        26*(-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        12*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2) + 12*(p1p3*p1p3) - 
                        17*(me*me)*(kp1 - kp3 + p1p3) - 
                        6*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                        16*(-kp1 + me*me + p1p2 - p1p3)*
                       (kp1 - kp3 + p1p3) + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                       (kp1 - kp3 + p1p3) + 
                        2*pow(kp1 - kp3 + p1p3,2) - 
                        2*p1p3*
                       (kp1 - kp3 + 12*(me*me) + 
                       4*(-kp1 - kp2 + p1p2) - 
                       11*(-kp1 + me*me + p1p2 - p1p3) + p1p3) + 
                        2*p1p2*
                        (5*(me*me) - 
                        11*(-kp1 + me*me + p1p2 - p1p3) - 4*p1p3 - 
                        3*(kp1 - kp3 + p1p3))) + 
                      pow(me,4)*
                       (-7*kp1 - 7*kp2 - 18*(me*me) - 18*p1p2 - 
                        3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        7*p1p3 + 
                        2*(-15*kp1 - kp2 + 6*(kp1 + kp2 - kp3) + 
                        me*me + 17*p1p2 + 
                        4*(-kp1 + me*me + p1p2 - p1p3) + 
                        2*(kp1 - kp3 + p1p3)))))/
                  (2*(me*me) - 2*(kp1 - kp3 + p1p3)) - 
                 (512*(8*pow(me,6) + 
                      me*me*(-kp1 - kp2 + p1p2)*
                       (-kp1 + me*me + p1p2 - p1p3) - 
                      2*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                      2*(-kp1 - kp2 + p1p2)*
                       (-kp1 + me*me + p1p2 - p1p3)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                      4*(-kp1 + me*me + p1p2 - p1p3)*
                       pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) \
- me*me*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                      4*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                       (kp1 - kp3 + p1p3) - 
                      4*(-kp1 - kp2 + p1p2)*
                       (-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) \
+ me*me*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                      2*(p1p2*p1p2)*(-kp1 - kp2 + p1p2 + p1p3) + 
                      pow(me,4)*
                       (kp2 - kp3 + 2*p1p2 + 2*(-kp1 - kp2 + p1p2) + 
                        2*(-kp1 + 2*(me*me) + 5*p1p2 + 
                        3*(-kp1 - kp2 + p1p2) - 2*p1p3) - p1p3 - 
                        3*(kp1 - kp3 + p1p3)) + 
                      p1p2*(4*(me*me)*(-kp1 - kp2 + p1p2) + 
                        2*pow(-kp1 - kp2 + p1p2,2) + 
                        (3*(me*me) + 2*(-kp1 - kp2 + p1p2) - 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        (-kp1 + me*me + p1p2 - p1p3) - 
                        2*(-kp1 - kp2 + p1p2)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        me*me*(kp1 - kp3 + p1p3) + 
                        p1p3*
                        (6*(me*me) - 6*(-kp1 + me*me + p1p2 - p1p3) - 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        6*(kp1 - kp3 + p1p3))) + 
                      me*me*(2*(p1p2*p1p2) + 
                        4*(me*me)*(-kp1 - kp2 + p1p2) + 
                        2*pow(-kp1 - kp2 + p1p2,2) + 
                        3*(me*me)*(-kp1 + me*me + p1p2 - p1p3) + 
                        3*(-kp1 - kp2 + p1p2)*
                        (-kp1 + me*me + p1p2 - p1p3) - 
                        me*me*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*(-kp1 - kp2 + p1p2)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        6*(-kp1 + me*me + p1p2 - p1p3)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2) + 
                        2*p1p2*
                        (kp2 - kp3 + me*me + 5*(-kp1 - kp2 + p1p2) - 
                        p1p3) + me*me*(kp1 - kp3 + p1p3) - 
                        (-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                        4*(-kp1 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) + 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) - 
                        p1p2*(-kp3 + me*me + p1p2 + 2*p1p3) + 
                        p1p3*
                        (3*(me*me) - 
                        8*(-kp1 + me*me + p1p2 - p1p3) - 
                        8*(kp1 - kp3 + p1p3)) - 
                        p1p3*
                        (3*(-kp1 + me*me + p1p2 - p1p3) - 
                        2*(kp1 - kp3 + p1p3))) + 
                      p1p3*(2*(kp1 - me*me)*(kp1 - kp3 + p1p3) + 
                        (-kp1 + me*me + p1p2 - p1p3)*
                        (me*me + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        6*(kp1 - kp3 + p1p3)))))/(me*me + p1p2))/(2.*kp2) \
- ((512*(8*pow(me,6) - me*me*(-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3) - 
                     4*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     2*(-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     2*pow(-kp1 + me*me + p1p2 - p1p3,2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     2*(p1p2*p1p2)*
                      (-2*kp1 - 2*kp2 + kp3 + me*me + 2*p1p2 - p1p3) + 
                     me*me*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                     4*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) - 
                     4*(-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) \
+ me*me*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) + 
                     6*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) + 
                     4*(p1p3*p1p3)*(kp1 - kp3 + p1p3) + 
                     p1p3*((me*me + 
                       2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        (-kp1 + me*me + p1p2 - p1p3) - 
                        2*(-kp3 + p1p3)*(kp1 - kp3 + p1p3)) + 
                     pow(me,4)*
                      (2*kp1 + kp2 - 2*kp3 - me*me + p1p2 + 
                        2*(-kp1 - kp2 + p1p2) - 
                        3*(-kp1 + me*me + p1p2 - p1p3) + p1p3 + 
                        2*(2*kp1 + kp2 - 2*kp3 + 3*p1p2 + 
                        3*(-kp1 - kp2 + p1p2) + 2*p1p3)) + 
                     p1p2*(4*(me*me)*(-kp1 - kp2 + p1p2) + 
                        2*pow(-kp1 - kp2 + p1p2,2) + 
                        (me*me - 
                       6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        (-kp1 + me*me + p1p2 - p1p3) - 
                        2*(-kp1 - kp2 + p1p2)*p1p3 + 
                        3*(me*me)*(kp1 - kp3 + p1p3) + 
                        2*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                        2*p1p3*(kp1 - kp3 + p1p3) + 
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (6*(me*me) - 2*p1p3 - 6*(kp1 - kp3 + p1p3))) + 
                     me*me*(2*(p1p2*p1p2) + 
                        4*(me*me)*(-kp1 - kp2 + p1p2) + 
                        2*pow(-kp1 - kp2 + p1p2,2) - 
                        p1p2*
                        (-kp3 + me*me + p1p2 + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)) + 
                        me*me*(-kp1 + me*me + p1p2 - p1p3) - 
                        (-kp1 - kp2 + p1p2)*
                        (-kp1 + me*me + p1p2 - p1p3) + 
                        3*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        6*(-kp1 + me*me + p1p2 - p1p3)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*(p1p3*p1p3) + 3*(me*me)*(kp1 - kp3 + p1p3) + 
                        3*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                        4*(-kp1 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) - 
                        11*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) + 
                        2*p1p2*
                        (2*kp1 + kp2 - 2*kp3 - p1p2 + 
                        5*(-kp1 - kp2 + p1p2) + p1p3) - 
                        p1p3*
                        (-kp1 + me*me + p1p2 - p1p3 + 
                        2*(kp1 - kp3 + p1p3)) - 
                        p1p3*
                        (me*me + 2*(-kp1 - kp2 + p1p2) - 
                        2*(-kp1 + me*me + p1p2 - p1p3) + 
                        4*(kp1 - kp3 + p1p3)))))/(me*me + p1p2) + 
                (256*(20*pow(me,6) + 
                     4*(-(kp1*kp2*(me*me)) - kp1*kp3*(me*me) + 
                        kp2*pow(me,4) + kp3*pow(me,4) + 
                        p1p2*
                        (-(kp2*(2*kp3 + me*me)) + 
                        kp3*(2*kp3 + 3*(me*me)) + 
                        2*(-kp2 + me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        p1p3,2)) - 
                        2*(kp1*kp1)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*kp1*kp2*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*kp1*kp3*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        4*kp1*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*kp3*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*(p1p2*p1p2)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*(me*me)*
                        pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2) + 
                        (2*(kp2*kp2) + kp3*(me*me) - 
                        kp2*(2*kp3 + 3*(me*me)) + 
                        2*(-2*kp1 + kp2 + me*me)*
                       (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                       p1p3,2))*p1p3 - 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (p1p3*p1p3)) + 
                     pow(me,4)*
                      (18*kp2 - 18*(me*me) - 2*p1p2 - 
                        7*(-kp1 + me*me + p1p2 - p1p3) + 
                        11*(kp1 - kp3 + p1p3) + 
                        2*(kp1 + kp2 - kp3 + 4*p1p2 + 
                        4*(-kp1 - kp2 + p1p2) + 
                        5*(-kp1 + me*me + p1p2 - p1p3) + 2*p1p3 - 
                        7*(kp1 - kp3 + p1p3))) + 
                     me*me*(-14*pow(me,4) - 4*(p1p2*p1p2) - 
                        24*(me*me)*(-kp1 - kp2 + p1p2) - 
                        12*pow(-kp1 - kp2 + p1p2,2) + 
                        9*(me*me)*(-kp1 + me*me + p1p2 - p1p3) + 
                        18*(-kp1 - kp2 + p1p2)*
                        (-kp1 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 + me*me + p1p2 - p1p3,2) + 
                        12*(me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        12*(-kp1 - kp2 + p1p2)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        14*(-kp1 + me*me + p1p2 - p1p3)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        16*(p1p3*p1p3) + 
                        45*(me*me)*(kp1 - kp3 + p1p3) + 
                        22*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                        8*(-kp1 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) - 
                        22*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) - 
                        18*pow(kp1 - kp3 + p1p3,2) + 
                        2*p1p2*
                        (-kp1 + 4*(me*me) + p1p2 + 
                        12*(-kp1 - kp2 + p1p2) - 
                        6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        p1p3 + 5*(kp1 - kp3 + p1p3)) - 
                        p1p3*
                        (-34*(me*me) - 24*(-kp1 - kp2 + p1p2) + 
                        2*(-kp1 + me*me + p1p2 - p1p3) + 
                        34*(kp1 - kp3 + p1p3)) - 
                        4*((-kp1 - kp2 + p1p2)*
                        (-kp1 + me*me + p1p2 - p1p3) - 
                         2*(-kp1 + me*me + p1p2 - p1p3)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                         p1p2*
                        (-kp1 + me*me + p1p2 + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        p1p3) + (-kp1 + me*me + p1p2 - p1p3)*p1p3 + 
                         (-kp1 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) + 
                         2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                         (kp1 - kp3 + p1p3)))))/
                 (2*(me*me) - 2*(kp1 - kp3 + p1p3)))/(2.*kp1))/
           (2.*(-kp1 - kp2 + me*me + p1p2)) - 
          ((128*(12*pow(me,6) + 
                  4*(-(kp1*kp2*(me*me)) - kp1*kp3*(me*me) + 
                     kp2*pow(me,4) + kp3*pow(me,4) + 
                     p1p2*(-(kp2*(2*kp3 + me*me)) + 
                        kp3*(2*kp3 + 3*(me*me)) + 
                        2*(-kp2 + me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                        2)) - 
                     2*(kp1*kp1)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     2*kp1*kp2*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     2*kp1*kp3*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     4*kp1*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     2*kp3*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     2*(p1p2*p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     2*(me*me)*
                      pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                     (2*(kp2*kp2) + kp3*(me*me) - 
                        kp2*(2*kp3 + 3*(me*me)) + 
                        2*(-2*kp1 + kp2 + me*me)*
                        (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - 
                        p1p3,2))*p1p3 - 
                     2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (p1p3*p1p3)) + 
                  pow(me,4)*
                   (9*(kp1 + kp2 - kp3) + 7*kp3 + 8*(me*me) + 
                     24*(-kp1 - kp2 + p1p2) - 
                     16*(-kp1 + me*me + p1p2 - p1p3) + 4*p1p3 + 
                     2*(-15*(kp1 + kp2 - kp3) + 6*kp3 - 3*(me*me) + 
                        5*(-kp1 - kp2 + p1p2) + 
                        4*(-kp1 + me*me + p1p2 - p1p3) + 4*p1p3)) + 
                  me*me*(-4*(kp1*kp1) - 4*kp1*(kp1 + kp2 - kp3) + 
                     12*pow(kp1 + kp2 - kp3,2) + 26*kp1*kp3 - 
                     14*(kp1 + kp2 - kp3)*kp3 - 6*(kp3*kp3) + 
                     20*kp1*(me*me) + 17*(kp1 + kp2 - kp3)*(me*me) - 
                     25*kp3*(me*me) - 6*pow(me,4) + 
                     36*kp1*(-kp1 - kp2 + p1p2) + 
                     4*(kp1 + kp2 - kp3)*(-kp1 - kp2 + p1p2) - 
                     14*kp3*(-kp1 - kp2 + p1p2) - 
                     34*(me*me)*(-kp1 - kp2 + p1p2) - 
                     8*pow(-kp1 - kp2 + p1p2,2) + 
                     2*(-8*kp1 + 9*kp3 + 6*(me*me) + 
                        18*(-kp1 - kp2 + p1p2))*
                      (-kp1 + me*me + p1p2 - p1p3) - 
                     12*pow(-kp1 + me*me + p1p2 - p1p3,2) - 
                     2*(6*kp1 + 10*(kp1 + kp2 - kp3) - 7*kp3 + 
                        10*(me*me) - 4*(-kp1 - kp2 + p1p2) - 
                        18*(-kp1 + me*me + p1p2 - p1p3))*p1p3 - 
                     4*((-kp1 - kp2 + p1p2)*
                        (-kp1 + me*me + p1p2 - p1p3) - 
                        2*(-kp1 + me*me + p1p2 - p1p3)*
                         (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        p1p2*
                         (-kp1 + me*me + p1p2 + 
                         2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                         p1p3) + (-kp1 + me*me + p1p2 - p1p3)*p1p3 + 
                        (-kp1 + me*me + p1p2 - p1p3)*
                        (kp1 - kp3 + p1p3) + 
                        2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                         (kp1 - kp3 + p1p3)))))/
              ((me*me + p1p2)*(2*(me*me) - 2*(kp1 - kp3 + p1p3))) - 
             (128*(-64*pow(me,6) + 
                  pow(me,4)*
                   (-19*kp1 - 25*kp2 + 16*kp3 + 28*(me*me) + 38*p1p2 - 
                     36*p1p3 - 
                     2*(19*(me*me) + 11*p1p2 + 28*(-kp1 - kp2 + p1p2) - 
                        38*(-kp1 + me*me + p1p2 - p1p3) - 
                        23*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        23*p1p3 - 46*(kp1 - kp3 + p1p3))) + 
                  4*(6*(me*me)*(-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3) - 
                     4*(me*me)*(-kp1 - kp2 + p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     4*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     5*(-kp1 - kp2 + p1p2)*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     pow(-kp1 + me*me + p1p2 - p1p3,2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     3*(-kp1 + me*me + p1p2 - p1p3)*
                      pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) - 
                     4*(me*me)*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) + 
                     2*(me*me)*(-kp1 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) - 
                     4*(me*me)*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) - 
                     (-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) - 
                     7*(p1p3*p1p3)*(kp1 - kp3 + p1p3) + 
                     2*(p1p2*p1p2)*
                      (kp1 + kp2 - kp3 - me*me - p1p2 - 
                        2*(-kp1 - kp2 + p1p2) + 2*p1p3) - 
                     p1p3*(-4*(me*me)*(-kp1 - kp2 + p1p2) + 
                        (-3*(-kp1 - kp2 + 4*(me*me) + p1p2) + 
                        7*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                        (kp1 - kp3 + p1p3) + 
                        7*pow(kp1 - kp3 + p1p3,2) + 
                        (-kp1 + me*me + p1p2 - p1p3)*
                         (-kp1 - kp2 + kp3 + 3*(me*me) + p1p2 - p1p3 + 
                         3*(kp1 - kp3 + p1p3))) + 
                     p1p2*(-2*(me*me)*(-kp1 - kp2 + p1p2) + 
                        (-3*(me*me) - 4*(-kp1 - kp2 + p1p2) + 
                        11*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3))*
                         (-kp1 + me*me + p1p2 - p1p3) - 
                        19*(me*me)*
                         (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                        10*(-kp1 - kp2 + p1p2)*
                         (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        2*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,
                         2) - 2*(p1p3*p1p3) + 
                        3*(me*me)*(kp1 - kp3 + p1p3) + 
                        6*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                         (kp1 - kp3 + p1p3) + 
                        p1p3*(5*(me*me) + 2*(-kp1 - kp2 + p1p2) - 
                         2*(-kp1 + me*me + p1p2 - p1p3) + 
                         4*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                         5*(kp1 - kp3 + p1p3)))) + 
                  me*me*(-10*pow(me,4) + 14*(p1p2*p1p2) - 
                     128*(me*me)*(-kp1 - kp2 + p1p2) - 
                     56*pow(-kp1 - kp2 + p1p2,2) - 
                     31*(me*me)*(-kp1 + me*me + p1p2 - p1p3) + 
                     12*(-kp1 - kp2 + p1p2)*
                      (-kp1 + me*me + p1p2 - p1p3) + 
                     12*pow(-kp1 + me*me + p1p2 - p1p3,2) + 
                     85*(me*me)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                     106*(-kp1 - kp2 + p1p2)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     54*(-kp1 + me*me + p1p2 - p1p3)*
                      (-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) - 
                     50*pow(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3,2) + 
                     20*(p1p3*p1p3) + 7*(me*me)*(kp1 - kp3 + p1p3) + 
                     98*(-kp1 - kp2 + p1p2)*(kp1 - kp3 + p1p3) - 
                     30*(-kp1 + me*me + p1p2 - p1p3)*(kp1 - kp3 + p1p3) - 
                     60*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                      (kp1 - kp3 + p1p3) - 26*pow(kp1 - kp3 + p1p3,2) - 
                     2*p1p2*(8*(me*me) + 29*(-kp1 - kp2 + p1p2) + 
                        13*(-kp1 + me*me + p1p2 - p1p3) - 
                        18*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        9*p1p3 + 2*(kp1 - kp3 + p1p3)) + 
                     p1p3*(15*(me*me) + 68*(-kp1 - kp2 + p1p2) - 
                        46*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                        58*(kp1 - kp3 + p1p3)) - 
                     4*(p1p2*(2*(-kp1 - kp2 + p1p2) - 
                         3*(-kp1 + me*me + p1p2 - p1p3) - 
                         7*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + 
                         p1p3 + 3*(kp1 - kp3 + p1p3)) - 
                        2*((-kp3 + me*me + p1p2 - 2*(-kp1 - kp2 + p1p2))*
                         p1p3 + 
                         2*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3)*
                         (-kp2 - kp3 + p1p2 + p1p3) + 
                         (-kp1 + me*me + p1p2 - p1p3)*
                         (-kp2 - kp3 + p1p2 - 
                         3*(-kp1 - kp2 + kp3 + me*me + p1p2 - p1p3) + p1p3\
))))))/pow(2*(me*me) - 2*(kp1 - kp3 + p1p3),2))/(2.*kp1))/(2.*kp3)))/64.;
}
