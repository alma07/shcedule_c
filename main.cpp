#include <QCoreApplication>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <conio.h>

using namespace std;

/*#include <middle_begin.cpp>*/

const int N = 3;
const int V = 10;
const int Z = 5;
const int L = 3;
const float epsilon = 0.3; //точность
/*константы*/
int iteration=1;        //для отметки логов (т.к. их несколько за одно выполнение программы)
char route_log[22]="C:\\output\\log\\middle_";
const int II=N;          //количество типов
//const int L=2;          //количество приборов
const int QQ=15;         //максимальное количество требований в партии
const char* route_til= "C:\\output\\time\\til.txt";
const char* route_tikl="C:\\output\\time\\tikl.txt";
const char* route_tiil="C:\\output\\time\\tiil.txt";
const char* route_tjq0l="C:\\output\\time\\tjq0l.txt";
const char* route_tjinl="C:\\output\\time\\tjinl.txt";
const char* route_middle="C:\\output\\middle.txt";
const char* route_P="C:\\output\\P.txt";
const char* route_R="C:\\output\\R.txt";

const int razmer=10;

/*матрицы времён*/
float t_i_l[II+1][L+1],          //обработки требования i-го типа на l-ом приборе
      t_ik_l[II+1][II+1][L+1],    //переналадки l-го приб. с обр-ки i-го типа на k-ый
      t_ii_l[II+1][L+1];         //первонач. наладки l-го приб. на обраб i-го типа

struct desision {
    int M[V];
    int A[V][V];

};
struct party { //партия
    int type;
    int m;
    int A[V];

};
//**********************переменные для верхнего уровня******************************
desision current,//текущее решение
        global,//глобально эффективное
        local,//эффективное для i типа
        star;//эффективное для заданного кол-ва партий
int s = 0, //шаг
    M = 100.0,
   // N = 3, //кол-во типов партий (надо читать динамически с формы)
    h = 0, //h'
    H = 0,
    j = 0,//i'
    cicle = 0; //переход по пунктам алгоритма
int I[N+1][2]; // 0 - кол-во требований данного типа
             // 1 - "0" - не было анализа, "1" - был
int m[N+1];//количество партий каждого типа
float fG[V+1];
float fL[V];
float fStar[V] ;
float f[V];
float delta_minus = 0; //дискретные градиенты
float delta_plus = 0;

int ii = 0, jj = 0, k = 0; //counter
int As[V][V]; //current A(S)
int As1[V][V]; //A(S+1)

bool right = true;
bool finish = false;
float t_z = 14; //время обработки партии (задано)
//****************************************

//**********************переменные для среднего уровня******************************
int z, //индекс текущей партии
    kz[Z],
    kq,
    k_mb,
    i_i; //номер типа требования
int m_i[V]; //mi' количество партий требований i  типа, размещенных в группах N_z  и в множестве Q
int h_mb;//номер партиии
int mi[V];
party N_z[Z][Z], Q[Z]; //множество партий
int NZ[Z][2],
    NZT[Z][2]; //множество номеров (идентификаторов) групп партий
bool correct_19 = false;
bool correct_20 = false;
bool finish_mb = false;
int Zz; //количество групп сформированных на начальном решении
//----------------эффективное решение-------------------------
int z_e, //индекс текущей группы партий, которая оптимизируется z'
    I_z[N+1],//(I)z
    I_ze[N+1],//(I)z'
    //номера партий
    Jp[3][3], //будут использоваться с индексами 1 и 2
    Ip[3][3], //будут использоваться с индексами 1 и 2
    P[3], //номера варианто решения
    l1, //сохранненное решение
    i_p[3],//i'p
    ip[3],//ip
    r_p[3],//r'p
    rp[3],//rp
    k_zp[3],//k'p
    kzp[3],//kp
    i_star,//i'*
    istar,//i*
    r_star,//r'*
    rstar,//r*
    h_e, //h'
    i_ud[3], i_dob[3],
    r_ud[3], r_dob[3];
int f_count = 0; //для заполнения функций
int kk = 0, hh = 0;//для незавершенного производства


float f2_z[V], f2_ze[V], //целевые функции
    leftGrad_z[3], leftGrad_ze[3],
    rightGrad_z[3], rightGrad_ze[3],
    G1, G2,
    fz_2p[3], fz2p[3];

party Nz_local[Z], Nz_e_local[Z]; //локально эффективное решение
party N_e[3], N_num[3];
party cur_N_ze[3][Z], cur_N_z[3][Z];

party N_z_effective[Z][Z]; //глобально эффективное решение

bool finish_em = false;
// нулевые столбцы и строки рассматривать не будут
int main_P[N+1][V+1], main_R[N+1][V+1]; //уточнить размерности
float tjq0l[L+1][V+1][QQ+1], til[V+1][V+1];

// нулевые столбцы и строки рассматривать не будут
int main_P_e[N+1][V+1], main_R_e[N+1][V+1]; //уточнить размерности
float tjq0l_e[L+1][V+1][QQ+1];

//------------------------------------------------------------
//****************************************

//***************************функции верхнего уровня*****************
int beginDesision(int t);
bool checkCorrect(int t);
int step_3(int t);
int step_4();
int step_7();
int step_8();
int step_10();
int step_11();
int step_12();
int step_14();
int step_15();
int step_16();
int step_17();
int step_18_19_20();
int step_22(int t);
int step_23_24(int t);
int step_25();
int step_26_27_28();
void rememberDesision();
void fixDesision(int value);
int find_f1();

//**************************************************

//***************************функции среднего уровня*****************
void begin_middle(desision check);
int mb_step_3(int);
int mb_step_5(int t, desision x);
int mb_step_6(party x, desision y, int t);
int mb_step_8(int t);
void  mb_step_7();
int mb_step_9();
int mb_step_10(int t, desision x);
int mb_step_11(int t);
int mb_step_12(int t, desision x);
int mb_step_13(int t, party x);
int mb_step_14(int t);
int mb_step_15_16(int t);
int mb_step_17_18(int t, desision x);
int mb_step_19(int t);
int mb_step_20();
int mb_step_21_22(int t, desision x);
int mb_step_23_24(int t,desision x);
bool check_condition(int condition);
bool check_condition_e(int);
bool check_condition_num(int);
float find_f2(int nzp, int nnzp, int *nj, int nz);
float find_f2_e(int nzp, int nnzp, int *nj, int nz);
void PartyCount();
void ReWatchGroup();
int partyInGroupCount(int x);
void writeGroupLog();
void writeEffGroupLog();

//----------------эффективное решение-------------------------
void effective_middle();
void readFromFile(int param, int, int, int); // 0 - P, 1 - R, 2 - tjq0l, 3 - til
//вернет количество строк в файле
int writeToFile(const char*, int x, int m); //m = 0(N_z), m=1(cur_N_z), m=2(cur_N_ze)
void WriteFunc2ToFile();
void WriteGroupSwitchLog(int, int, int, int, int, int, char*);
int find_j11(); //1 - в степени
int find_j12(); //1 - в степени
int find_j21(int t); //2 - в степени
int find_j22(int t); //2 - в степени
int e_step_2();
int e_step_3();
int e_step_5(int, int);
void form_I_z(int one, int two);
int e_step_6(int t);
int e_step_7(int t);
int e_step_8();
int e_step_10();
int e_step_11();
int e_step_12();
int e_step_13(int, int);
int e_step_14(int);
int e_step_15();
int e_step_16();
int e_step_18();
int e_step_19(int num);
int e_step_20(int);
int e_step_21_22(int);
int e_step_23(int);
int e_step_24_25_26(int);
int e_step_27();
int e_step_37();
int e_step_38();
int e_step_39();
int e_step_40();
int e_step_41_42_43(int t);
int e_step_44();
int e_step_45();
int e_step_46_47();
int e_step_48_49();
int e_step_50(int);
int e_step_51();
int e_step_52();
int e_step_53(int);
int e_step_54_55(int num);
int e_step_56();
int e_step_58();
int e_step_63();
int e_step_64_65();
int e_step_66_67_68();
int e_step_69();
int e_step_70();
int e_step_71();
int e_step_72_73();

//------------------------------------------------------------
//**************************************************

//*********************Функции нижнего уровня**********
void start_part(int **P, int **R, int n_z, int n_p_z);
float objective(int **P, int **R, int n_z, int n_p_z);
void start_time (int **P, int **R, int n_z, int n_p_z);
int low_level(int x); //количество строк во входном файле
//*****************************************************
int main(int argc, char *argv[])
{
    ofstream out_file;
    ofstream out_q;
    ofstream out_eff_file;
    ofstream out_log_eff_file;
    out_log_eff_file.open("C://output//effective_log.txt", ios::app);
    out_file.open("C:\\output\\group.txt", ios::out);
    out_q.open("C:\\output\\q.txt", ios::out);
    out_eff_file.open("C:\\output\\eff_group.txt", ios::out);
    out_file.close();
    out_q.close();
    out_eff_file.close();
    out_log_eff_file.close();
    //начальная инициализация матриц с нижнего уровня
    for (s=0;s<=N;s++){
        for (ii=0;ii<=Z;ii++){
            main_P[s][ii] = 0;
            main_R[s][ii] = 0;
            main_P_e[s][ii] = 0;
            main_R_e[s][ii] = 0;
        }
    }

    for (s=0;s<=V;s++){
        for (ii=0;ii<=V;ii++){
            til[s][ii] = 0;
            til[s][ii] = 0;
        }
    }

    for (s=0;s<=L;s++){
        for (ii=0;ii<=V;ii++){
            for (j=0;j<=QQ;j++){
                tjq0l[s][ii][j] = -1;
                tjq0l_e[s][ii][j] = -1;
            }
        }
    }

   //время обработки i  типа на l приборе (задано)
    til[1][1] = 1;  til[1][2] = 2;  til[1][3] = 1;
    til[2][1] = 1;  til[2][2] = 2;  til[2][3] = 1;
    til[3][1] = 1;  til[3][2] = 2;  til[3][3] = 1;
    til[4][1] = 1;  til[4][2] = 2;  til[4][3] = 1;

//----------------------------------------------
//задание количества партий
    /*for (ii=1;ii<N+1;ii++){
        I[ii][0] = N*ii-1;
    }*/
    I[1][0] = 7; I[2][0] = 7;
    I[3][0] = 7; //I[4][0] = 5;
    //I[5][0] = 5; I[6][0] = 6;
//--------------------------------
    for(ii=0;ii<V;ii++){
        for(s=0;s<V;s++){
            As[ii][s] = 0;
        }
    }

    //step 1
    s = 0;
    j = 1; I[1][1] = 1;//i'=1;
    m[1] = 2;//mi' = 2;
    for(ii=2; ii<N+1; ii++){
        m[ii] = 2;
    }
    fG[s] = M;

    //step 2
    for(ii=2; ii<N+1; ii++){
        beginDesision(ii);
        checkCorrect(ii);
    }

    int cicle = 3;
    while (true){ //основная программа
        switch(cicle){
            case 3: {
                cicle = step_3(j);//формирование начального состава партий для i'
                rememberDesision(); //формирование текущего решения
                if (cicle == 26){
                    fixDesision(0); //инициализация глоабльного решния как локального
                }
                break;}
            case 4: {
                 cicle = step_4();//фиксация текущего решения как глобального
                 break; }
            case 5: //передача текущего решения на 2 уровень иерархии
                begin_middle(current);
                cicle = 6;
                break;
            case 6:
                //получение решения с верхнего уровня
                cicle = 7;
                //f[s] = rand()%100;
                f[s] = find_f1();
                s++;
                break;
            case 7:
                //фиксация текущего как локального эффективного
                cicle = step_7();
                break;
            case 8:
                //фиксация текущего как локального эффективного
                cicle = step_8();
                break;
            case 9:
                h = 2;
                cicle = 10;
                break;
            case 10:
                cicle = step_10();
                break;
            case 11:
                cicle = step_11();
                break;
            case 12:
                cicle = step_12();
                break;
            case 13:
                //получение эффективного решения с другого уровня
                //f[s+1] = rand();
                f[s] = find_f1();
                //f[s] = rand()%50;
                cicle = 14;
                break;
            case 14: //анализ эффективности
                cicle = step_14();
                break;
            case 15:
                cicle = step_15();
                break;
            case 16:
                cicle = step_16();
                break;
            case 17:
                cicle = step_17();
                break;
            case 18:
                cicle = step_18_19_20();
                break; //step 18-19-20
            //case 19: break;
            //case 20: break;
            case 21:
                cicle = 22;
                m[j]++;
                break;
            case 22:
                cicle = step_22(j);
                break;
            case 23:
                cicle = step_23_24(j);
                break;
            //case 24: break;
            case 25:
                cicle = step_25();
                break;
            case 26:
                cicle = step_26_27_28();
                break;
            //case 27: break;
            //case 28: break;
        }
        if (finish){//проверка окончания
            break;
        }
    }//конец основной программы
    std::cout<<'1';
    return 0;
}

//тексты функций
//------------------------------------------------------------
void WriteFunc2ToFile(){
    ofstream out_file;
    out_file.open("C:\\output\\func2.txt", ios::app);

    out_file << "f2_z[0]="<<f2_z[0] << "\t";
    out_file << "f2_ze[0]="<<f2_ze[0] << "\t";
    out_file << "fz_2p[1]="<<fz_2p[1] << "\t";
    out_file << "fz_2p[2]="<<fz_2p[2] << "\t";
    out_file << "fz2p[1]="<<fz2p[1] << "\t";
    out_file << "fz2p[2]="<<fz2p[2] << "\t";
    out_file << "\n";
}

//m = 0(N_z), m=1(cur_N_z), m=2(cur_N_ze)
int writeToFile(const char* path, int x, int m){ //x  - индекс группы
   //ofstream out_file;
   //out_file.open(path, "ios::out");

   ofstream out_file(path);
   int w = 0, ww = 0;
   int str = 0;
   if (m){
       if(m == 1){ //(cur_N_z)
           for (w = 0; w < Z; w++){
               if (cur_N_z[x][w].m != 0) {
                   str++;
                   out_file << cur_N_z[x][w].type << "\t";
                   out_file << cur_N_z[x][w].m << "\t";
                   for (ww = 1; ww < V; ww++){
                       if (cur_N_z[x][w].A[ww] != 0){
                           out_file << cur_N_z[x][w].A[ww] << "\t";
                       }
                   }
                   out_file << "\n";
               }
           }

       } else { //(cur_N_ze)
           for (w = 0; w < Z; w++){
               if (cur_N_ze[x][w].m != 0) {
                   str++;
                   out_file << cur_N_ze[x][w].type << "\t";
                   out_file << cur_N_ze[x][w].m << "\t";
                   for (ww = 1; ww < V; ww++){
                       if (cur_N_ze[x][w].A[ww] != 0){
                           out_file << cur_N_ze[x][w].A[ww] << "\t";
                       }
                   }
                   out_file << "\n";
               }
           }
       }

   } else { //N_z
       //ReWatchGroup();
        for (w = 0; w < Z; w++){
            if (N_z[x][w].m != 0) {
                str++;
                //cout << N_z[x][w].type;
                out_file << N_z[x][w].type << "\t";
                out_file << N_z[x][w].m << "\t";
                for (ww = 1; ww < V; ww++){
                    if (N_z[x][w].A[ww] != 0){
                        out_file << N_z[x][w].A[ww] << "\t";
                    }
                }
                out_file << "\n";
            }
        }
   }
   out_file.close();
   return str;
}


//принимает i - тип требования
int beginDesision(int t){ //начальное решение
    int q = 0;
    m[t] = 2;
    for (q=2; q<=m[t]; q++){
        As[t][q] = 2;
    }
    As[t][1] = I[t][0];
    for (q=2; q<=m[t]; q++){
        As[t][1] = As[t][1] - As[t][q];
    }
    return 0;
}

//принимает i - тип требования
bool checkCorrect(int t){ //проверка корректности решения
    int q = 0, w = 0;

    if(As[t][1] < 2){//формируется 1 партия
        m[t] = 1;
        for (q=0; q<V; q++){
            As[t][q] = 0;//обнуляем всю строку
        }
        As[t][1] = I[t][0];
        I[t][1] = 1;
        return false; //решение не корректное
    }
    for (q=2; q<=m[t]; q++){
        if(As[t][1] < As[t][q]){
            I[t][1] = 1;
            m[t] = 1;
            for (w=0; w<V; w++){
                As[t][w] = 0;//обнуляем всю строку
            }
            As[t][1] = I[t][0];
            return false; //нет необходимости просматривать дальше
        }
    }
    return true; //решение корректное
}

int step_3(int t){
    beginDesision(t);
    if (checkCorrect(t)){
        return 4;
    } else {        
        return 26;
    }
}

int step_4(){
    fixDesision(0);
    return 5;
}

int step_7(){
    fixDesision(1);
    int q = 0;
    for(q=0;q<V;q++){
        fL[q] = f[q];
    }
    return 8;
}

int step_8(){
    fixDesision(2);
    int q = 0;
    for(q=0;q<V;q++){
        fStar[q] = f[q];
    }
    return 9;
}

int step_10(){
    int q = 0;
    int w = 0;
    for (w = 0; w<V; w++){
        if (w != j){
            for (q = 0; q<V; q++){
               As1[w][q] = As[w][q];
            }
        }
    }
    if (As1[j][h] == 0){
        As1[j][h] = As[j][h]+1;
    } else {
        As1[j][h] = As1[j][h]+1;
    }
    for (q=2;q<=m[j];q++){
        if(q != h){
             As1[j][q] = As[j][q];
        }
    }
    //As1[j][1] = m[j];
    As1[j][1] = I[j][0];
    for (q=2;q<=m[j];q++){
        As1[j][1] = As1[j][1] - As1[j][q];
    }
    return 11;
}

int step_11(){ //проверка выполнения условия окончания процесса изменения состава партий требований  i'-го типа
    if(As1[j][h] > As1[j][1]){
        return 18; //формирование прекращается
    } else {
        return 12;
    }
}

int step_12(){
    if(As1[j][h] <= As1[j][1]){
        desision next;
        int w = 0, ww = 0;
        for (w = 0; w<V; w++){
            next.M[w] = m[w];
            for (ww=0; ww<V; ww++){
                next.A[w][ww] = As1[w][ww];
            }
        }
        //решение As1 передается на верхний уровень
        begin_middle(next);
        return 13;
    }
}

int step_14(){
    delta_minus = delta_plus = f[s] - f[s-1];
    //s++;
    return 15;
}

int step_15(){
    int q = 0, w = 0;
    if (delta_minus < 0){
        //фиксация решения
        for(q=0;q<V;q++){
            fStar[q] = f[q];
            star.M[q] = m[q];
            for(w=0;w<V;w++){
                star.A[q][w] = As1[q][w];
            }
        }

        return 10;
    }
    return 16;
}

int step_16(){
    int q = 0, w = 0;
    if (delta_plus > 0){
        //решение не рассматривается
        for(q=0;q<V;q++){
            for(w=0;w<V;w++){
                As1[q][w] = 0;
            }
        }
    }
    return 17;
}

int step_17(){
    h++;
    if (h <= m[j]){
        return 10;
    } else {
    //процесс формирования составов партий требований i' -го типа в количестве mi'  должен быть завершен
        return 18;
    }
}

int step_18_19_20(){
    int q = 0, w = 0;
    if (fStar[s] <fL[s]){
        for(q=0;q<V;q++){
            local.M[q] = star.M[q];
            fL[q] = fStar[q];
            for(w=0;w<V;w++){
                local.A[q][w] = star.A[q][w];
            }
        }

        //fL = fStar;
        return 21;
    } else {
        //решение не фиксируется
        return 21;
    }
}


int step_22(int t){
    int q = 0;
    for (q=2; q<=m[t]; q++){
        As[t][q] = 2;
    }
    As[t][1] = I[t][0];
    for (q=2; q<=m[t]; q++){
        As[t][1] = As[t][1] - As[t][q];
    }
    //rememberDesision();
    return 23;
}

int step_23_24(int t){
    if (checkCorrect(t)){
        //решение As передается на верхний уровень
        rememberDesision();
        begin_middle(current);
        //получение эффективного решения        
        f[s] = find_f1();
        s++;
        return 8;
    } else {
        //исследование эффективности завершено
        return 25;
    }
}


int step_25(){
    int q = 0, w = 0;
    //if (fL[s+1] < fG[s]){
    if (fL[s] < fG[s-1]){
        for(q=0;q<V;q++){
            global.M[q] = local.M[q];
            for(w=0;w<V;w++){
                global.A[q][w] = local.A[q][w];
            }
        }
        for(q=0;q<V;q++){
            fG[q+1] = fL[q];
        }
    }
    return 26;
}

int step_26_27_28(){
    int q = 0, w = 0;
    bool find = false;
    for (q=1;q<N+1;q++){
        if (I[q][1] == 0){ //не было еще анализа
            find = true;
            j = q;
            m[j] = 2;
            I[q][1] = 1;
            break; //останавливаем цикл на первом
        }
    }

    if (find){ //хоть одно еще не проанализировано
        for(q=0;q<V;q++){
            current.M[q] = global.M[q];
            for(w=0;w<V;w++){
                current.A[q][w] = global.A[q][w];
            }
        }
        return 5;
    } else { //уже все пройдены
        //фиксируется глобальное решение
        finish = true;
    }
}

void rememberDesision(){
    int q = 0, w = 0;
    for(q=0;q<V;q++){
        current.M[q] = m[q];
        for(w=0;w<V;w++){
            current.A[q][w] = As[q][w];
        }
    }
}

void fixDesision(int value){ //0 - глобальное, 1 - локальное для i' типа,
    //2 - локальное для заданного количества партий
    int q = 0, w = 0;
    if(value){
        if (value == 1){
            for (q=0;q<V;q++){
                local.M[q] = current.M[q];
                for (w=0;w<V;w++){
                    local.A[q][w] = current.A[q][w];
                }
            }
        } else {
            for (q=0;q<V;q++){
                star.M[q] = current.M[q];
                for (w=0;w<V;w++){
                    star.A[q][w] = current.A[q][w];
                }
            }
        }
    } else {
        for (q=0;q<V;q++){
            global.M[q] = current.M[q];
            for (w=0;w<V;w++){
                global.A[q][w] = current.A[q][w];
            }
        }
    }

}

int find_f1(){ //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
int w = 0, ww = 0, www = 0;
int sum_1[2];
sum_1[0] = sum_1[1] = 0;
for (w = 0; w<=N; w++){
    sum_1[0] += I[w][0]; //количество всех партий всех типов
}

for (w = 0; w < Z; w++){
    for (ww = 0; ww < kz[w]; ww++){
        for(www = 1; www < V; www++){
            sum_1[1] += N_z_effective[w][ww].A[www]; //здесь надо будет поменять на эффективное решение
        }
    }
}
return sum_1[0] - sum_1[1];
}
//------------------------------------------------------------

//-----------------функции среднего уровня----------------------
void begin_middle(desision check){
    int q = 0;
    int step = 0;

    for (q=0; q<Z; q++){
      NZ[q][0] = NZT[q][0] = q;
      NZ[q][1] = NZT[q][1]= 0; //не просмотрено

    }
    //очистка групп
    for (q = 0; q<Z; q++){
        kz[q] = 0;
        Q[q].m = 0;
        Q[q].type = 0;
        for (step = 0; step<Z; step++){
            N_z[q][step].type = 0;
            N_z[q][step].m = 0;
            Q[q].A[step] = 0;
        }
    }
    //начальная инициализация
    i_i = 1;
    kq = 0;
    k_mb = 0;

    for (q=0; q<V; q++){
        mi[q] = check.M[q];
    }
    //пустные множества N_z и Q
    z = 0; //z = 1; NZT(s+1) = NZT(s)\{z}
    kz[z] = 0;
    NZT[0][1] = 1; //просмотрено
    step = 3;
    finish_mb = false;
    while (true){ //начало основного цикла
        switch(step){
            case 3:
                step = mb_step_3(i_i);
                break;
            case 4:
               // h_mb = 1;
                h_mb = 1; //рассматривается 1 партия
                step = 5;
                break;
            case 5:
                step = mb_step_5(i_i, check);
                break;
            case 6:
               // не понятно что происходит и для чего надо
                step = mb_step_6(N_z[z][k_mb], check, i_i);
               break;
            case 7:
                step = 8;
                //N_z передается на 3 уровень, возвращается расписание
                mb_step_7();
                break;
            case 8:
                step = mb_step_8(i_i);
                break;
            case 9:
                step = mb_step_9();
                break;
            case 10:
                step = mb_step_10(i_i, check);
                break;
            case 11:
                step = mb_step_11(i_i);
                break;
            case 12:
                step = mb_step_12(i_i, check);
                break;
            case 13:
                step = mb_step_13(i_i, N_z[z][k_mb]);
                break;
            case 14:
                step = mb_step_14(i_i);
                break;
            case 15:
                step = mb_step_15_16(i_i);
                break;
            //case 16: break;
            case 17:
                step = mb_step_17_18(i_i, check);
                break;
           // case 18: break;
            case 19:
                step = mb_step_19(i_i);
                break;
            case 20:
                step = mb_step_20();
                break;
            case 21:
                step = mb_step_21_22(i_i, check);
                break;
            //case 22: break;
            case 23:
                step = mb_step_23_24(i_i, check);
                break;
            //case 24: break;
        }
        if (finish_mb){
            break;
        }

    }//конец основного цикла

    ReWatchGroup(); //!!!!!!!!!!!!!!!!!!!!!
    writeGroupLog();

    bool find = false;
    for (step = 1; step <Z; step++){
        if (partyInGroupCount(step)){ //если есть еще хоть одна группа есть смысл перестановок
            find = true;
            break;
        }
    }

    //переписываем эффективное из начального
    for (step = 0; step < Z; step++){
        for (q = 0; q < Z; q++){
            N_z_effective[step][q].m = N_z[step][q].m;
            N_z_effective[step][q].type = N_z[step][q].type;
            for (i_i = 0; i_i < V; i_i++){
                N_z_effective[step][q].A[i_i] = N_z[step][q].A[i_i];
            }
        }
    }

    if (find){ //если есть еще хоть одна группа есть смысл перестановок
    //эффективное решение
        effective_middle();
    } else {
        writeEffGroupLog();
    }
}

void writeEffGroupLog(){
    ofstream out_file;
    int w = 0, ww = 0, www = 0;
    out_file.open("C:\\output\\eff_group.txt", ios::app);
    for (w = 0; w< Z; w++){
        for (ww = 0; ww <Z; ww++){
            if(N_z_effective[w][ww].m > 0){
                out_file << N_z_effective[w][ww].type << "\t" << N_z_effective[w][ww].m << "\t";
                for (www = 1; www <V; www++){
                    if (N_z_effective[w][ww].A[www] > 0){
                        out_file << N_z_effective[w][ww].A[www]<< "\t";;
                    } else {
                        out_file << "\n";
                        break;
                    }
                }
            } else {
                break;
            }
        }
        out_file <<";"<< "\n";
    }
    out_file <<"********************"<<"\n";
    out_file.close();
}

void writeGroupLog(){
    ofstream out_file;
    ofstream out_q;
    int w = 0, ww = 0, www = 0;
    out_file.open("C:\\output\\group.txt", ios::app);
    out_q.open("C:\\output\\q.txt", ios::app);
    for (w = 0; w< Z; w++){
        if (Q[w].m > 0){
            out_q << Q[w].type << "\t" << Q[w].m <<"\t";
            for (ww = 1; ww <V; ww++){
                if (Q[w].A[ww] > 0){
                    out_q << Q[w].A[ww] << "\t";
                }
            }
            out_q <<"\n";
        }
        for (ww = 0; ww <Z; ww++){
            if(N_z[w][ww].m > 0){
                out_file << N_z[w][ww].type << "\t" << N_z[w][ww].m << "\t";
                for (www = 1; www <V; www++){
                    if (N_z[w][ww].A[www] > 0){
                        out_file << N_z[w][ww].A[www]<< "\t";;
                    } else {
                        out_file << "\n";
                        break;
                    }
                }
            } else {
                break;
            }
        }
        out_file <<";"<< "\n";
    }
    out_q <<"********************"<<"\n";
    out_file <<"********************"<<"\n";
    out_file.close();
    out_q.close();
}

void ReWatchGroup(){
    int w = 0, ww = 0;
    int www = 0;
    int counter = 0;

    for(w = 0; w < Z; w++){
        for(ww = 0;ww < Z; ww++){
            if (N_z[w][ww].m > 0 && N_z[w][ww].A[0] == 0 && N_z[w][ww].A[1] == 0){
                //N_z[w][ww].m --;
                N_z[w][ww].A[0] = 1; //метка - уже просмотрели
                for (www = 2; www < V; www++){
                    if ( N_z[w][ww].A[www-1] > 0){
                        counter++;
                    }
                    N_z[w][ww].A[www-1] = N_z[w][ww].A[www];
                }
                N_z[w][ww].m = counter++;
                counter = 0;
            } else {
                if (N_z[w][ww].m == 0){
                    N_z[w][ww].type = 0;
                    for (www = 0; www < V; www++){
                        N_z[w][ww].A[www] = 0;
                    }
                }
            }
        }
    }
}

int mb_step_3(int t){
    m_i[t] = 0;
    return 4;
}

int mb_step_5(int t, desision x){
    int w = 0;
    N_z[z][k_mb].type = t; //тип требований
    //N_z[z][k_mb].m = x.M[t];
    N_z[z][k_mb].m = mi[t];
    for (w=0;w<V;w++){
        //N_z[z][k_mb].A[w] = x.A[t][w];
    }

    kz[z]++;
    //k_mb = kz[z];
    k_mb++;
    mi[t] = 0;
    return 6;
    //return 7;
}

int mb_step_6(party x, desision y, int t){ //t - тип требования
    mi[t]++;//!!!!
   //N_z[z][k_mb].type  = t;
   //N_z[z][k_mb].m = mi[t];
   //N_z[z][k_mb].A[N_z[z][k_mb].m] = y.A[t][h_mb];

   int w = 0, ww = 0;
 /*  for (w = 0; w<Z; w++){
       for (ww = 0; ww<Z; ww++){
           if (N_z[w][ww].type == t && N_z[w][ww].m !=0){//партия включена в множество
               N_z[w][ww].m++;
               N_z[w][ww].A[N_z[w][ww].m] = N_z[w][ww].A[h_mb];
               N_z[w][ww].A[h_mb] = 0;
           }
       }
   }*/

   for (w = 0; w<Z; w++){
       if (N_z[z][w].type == t) {
           N_z[z][w].m = mi[t];
           N_z[z][w].A[mi[t]] = y.A[t][h_mb];
           // обнулим А после h_mb
           for (ww = (h_mb+1); ww<V; ww++){
                N_z[z][w].A[ww] = 0;
           }
           break;
       }
   }

   //mi[t]--;
    return 7;
}

void  mb_step_7(){
    int nz=0, nzp=0, nnzp=0, nj[Z];
    int w = 0, ww = 0;
    ReWatchGroup();
    //число партий в группе N_z1
    for(ww=0; ww<Z; ww++){
        if (N_z[z][ww].m != 0){
            nzp += (N_z[z][ww].m); //решение строится с количеством партий m+1
        }
    }

    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (N_z[z][ww].m != 0){
            nz++;
        }
    }
    //количество требований входящих в последнюю партию
    for(ww=kz[z]; ww>=0; ww--){
        if (N_z[z][ww].m != 0){
            for(w=V; w>=0; w--){
                if (N_z[z][ww].A[w]> 0){
                    nnzp = N_z[z][ww].A[w];
                    break;
                }
            }
            if(nnzp > 0){
                break;
            }
        }
    }

    int col = 0;
    //передача на нижний уровень Nz
    col = writeToFile(route_middle,z,0);
    low_level(col);
    //вычисление критериев
    //из расписания Piz
    //сначала передали на нижний уровень Nz
    readFromFile(0,0,nz, nzp);
    readFromFile(1,0,nz, nzp);
    readFromFile(2,0,nz, nzp);
}

int mb_step_8(int t){
    //проверка условия 19
    correct_19 = check_condition(19);
  /*  int y = rand()%10;
    if (y%2){
        correct_19 = true;
    } else {
        correct_19 = false;
    }*/
    if (correct_19){
        m_i[t]++;
        h_mb++;
    }
    return 9;
}

int mb_step_9(){
    //проверка условия 20
    correct_20 = check_condition(20);
  /*  int y = rand()%10;
    if (y%2){
        correct_20 = true;
    } else {
        correct_20 = false;
    }*/
    return 10;
}

int mb_step_10(int t, desision x){
    if(correct_19 && !correct_20){
        /*if (m_i[t] < mi[t]){
            return 6;
        }
        if (m_i[t] == mi[t]){
            return 23;
        }*/
        if (m_i[t] < x.M[t]){
            return 6;
        }
        if (m_i[t] == x.M[t]){
            return 23;
        }
    }
    return 11;
}

int mb_step_11(int t){    
    if(correct_20){
        //партии не могут быть больше размещены в NZ
        NZ[z][1] = 1; //просмотрено
        //определим следующую z
       /* for (w = 0; w<Z; w++){
            if (NZ[w][1] == 0){
                NZ[w][1] = 1;
                z = w;
                break;
            }
        }*/

    }
    return 12;
}

int mb_step_12(int t, desision x){
    int w = 0;
    bool empty = true; // по умолчанию пустое

    if(correct_19 && correct_20){
        for (w = 0; w<Z; w++){
            if (NZ[w][1] == 0){
                empty = false;
                break;
            }
        }
        if (m_i[t] < x.M[t] && !empty){
        //if (m_i[t] < mi[t] && !empty){
            for (w = 0; w<Z; w++){
                if(NZT[w][1] == 0){ //не просмотрено
                    z = w;
                    k_mb = kz[z];
                    NZT[w][1] = 1; //просмотрено
                    break; //самое первое
                }
            }
            return 5;
        }
        if (m_i[t] == x.M[t]){
        //if (m_i[t] == mi[t]){
            return 23;
        }
    }
    return 13;
}

int mb_step_13(int t, party x){
    if(!correct_19){
        int w = 0;
        //партия исключается из группы
        for (w = 0; w<Z; w++){
            if (N_z[z][w].type == t){
                 N_z[z][w].A[mi[t]] = 0;
                 N_z[z][w].m--;
                 break;
            }
        }
        //N_z[z][k_mb].A[mi[t]] = 0;
        mi[t]--;
        k_mb--;
        //N_z[z][k_mb].m = 0; //!!!
        //N_z[z][k_mb].A[t] = 0;

    }
    return 14;
}

int mb_step_14(int t){
    if (mi[t] == 0){ //партии не включаются в группу
        kz[z]--;
        N_z[z][k_mb].m = 0;
        //k_mb--; //!!!!!
    }
    return 15;
}

int mb_step_15_16(int t){
    int w = 0;
    bool empt = true; //пустое по умолчанию
    for (w=0;w<Z;w++){
        if (NZT[w][1] == 0){ //не просмотрено
            empt = false;
            break;
        }
    }
    if (empt){
        //партии в группах размещены быть не могут
        //проверяет условие на принадлежность незавершенному производству
    } else {
        z = w;
        k_mb = kz[z];
        NZT[w][1] = 1;
        return 5;
    }
    return 17;
}

int mb_step_17_18(int t, desision x){
    //принадлежит ли сформированная партия Q
    int w = 0, ww = 0;
    bool foundA = false;
    bool foundT = false;
    for (w=0;w<Z;w++){
        if (Q[w].type == t){ //тот же тип
            foundT = true;
            if (Q[w].A[h_mb] == x.A[t][h_mb]) {//та же партия
                foundA = true;
                break;
            }
            break;
        }

    }

    if (foundT){ //нашли тип
        if (!foundA){ // такую партию еще не добавляли
            Q[w].A[h_mb] = x.A[t][h_mb];
            Q[w].m++;
        }
    } else {
        Q[kq].type = t;
        Q[kq].m = 1;
        Q[kq].A[h_mb] = x.A[t][h_mb];
        kq++;
    }
    return 19;
}

int mb_step_19(int t){
    m_i[t]++;
    h_mb++;
    return 20;
}

int mb_step_20(){
    int w = 0;
    for (w=0;w<Z;w++){
        if (NZ[w][1] == 0){ //не просмотрено
            z = w;
            k_mb = kz[z];
            NZ[w][1] = 1;
            break;
        }
    }
    for (w=0;w<Z;w++){
        NZT[w][0] = NZ[w][0];
        NZT[w][1] = NZ[w][1];
    }
    return 21;
}

int mb_step_21_22(int t, desision x){
    //принадлежит ли сформированная партия N_z
    int w = 0, ww = 0;
    bool found = false;
    for (w=0;w<Z;w++){
        if (N_z[z][w].type == t){ //тот же тип
            for (ww = 0; ww < V; ww++){
                if (N_z[z][w].A[ww] == x.A[t][h_mb]){
                    found = true;
                    break;
                }
            }
            if (found){
                break;
            }
        }
    }

    //if (m_i[t] < mi[t]){
    if (m_i[t] < x.M[t]){
        if (found){
            return 6;
        } else {
            return 5;
        }
    }
    return 23;
   }
int mb_step_23_24(int t , desision x){
    int w = 0;
    //if (m_i[t] == mi[t]){
    if (m_i[t] == x.M[t]){
        //следующий тип требований
        i_i++;
        if (i_i <= N){ //N - количество типов требований
          for (w=0;w<Z;w++){
              NZT[w][0] = NZ[w][0];
              NZT[w][1] = NZ[w][1];
          }
          for (w=0;w<Z;w++){
              if (NZT[w][1] == 0){ //не просмотрено
                  z = w;
                  k_mb = kz[z];
                  NZT[w][1] = 1;
                  break;
              }
          }
          return 3;
        } else {
            finish_mb = true; // начальное решение сформировано
        }
    }
}

bool check_condition(int condition){
    //исследуется текущая группа z
    int w = 0, ww = 0, nzp = 0, nnzp = 0;
    float sum_1 = 0, sum_2[L+1];
    //число партий в группе N_z
    for(ww=0; ww<Z; ww++){
        if (N_z[z][ww].m != 0){
            nzp += (N_z[z][ww].m); //решение строится с количеством партий m+1
        }
    }
    //количество требований входящих в последнюю партию
    for(ww=kz[z]; ww>=0; ww--){
        if (N_z[z][ww].m != 0){
            for(w=V; w>=0; w--){
                if (N_z[z][ww].A[w]> 0){
                    nnzp = N_z[z][ww].A[w];
                    break;
                }
            }
            if(nnzp > 0){
                break;
            }
        }
    }
    for (ww=1; ww<=L; ww++){
        for (w=1; w<=kz[z];w++){
            sum_1 += til[w][ww]*main_P[w][nzp];
        }
        sum_2[ww] = tjq0l[ww][nzp][nnzp] + sum_1;
        sum_1 = 0;
    }

    sum_1 = sum_2[1]; //поиск максимума
    for (ww=2; ww<=L; ww++){
        if (sum_2[ww] > sum_1){
            sum_1 = sum_2[ww];
        }
    }
    if (condition == 20){
       // if (fabs(t_z - sum_1) == 0)
        if (fabs(t_z - sum_1) <= 0.001)
            return true;
        else return false;
    }
    if (condition == 19){
        //if (t_z >= sum_1)
        if(fabs(t_z - sum_1) <= 0.001 || (t_z - sum_1) > 0)
            return true;
        else return false;
    }
}

//условие 19 для эффективной части алгоритма для z'
bool check_condition_e(int t){ //номер 1 || 2
    //исследуется текущая группа z
    int w = 0, ww = 0, nzp = 0, nnzp = 0;
    int nz = 0;//число типов требований в текущей группе
    float sum_1 = 0, sum_2[L+1];
    //число партий в группе N_z
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[t][ww].m != 0){
            nzp += (cur_N_ze[t][ww].m); //решение строится с количеством партий m+1
            nz++;
        }
    }
    //количество требований входящих в последнюю партию
    for(ww=kz[t]; ww>=0; ww--){
        if (cur_N_ze[t][ww].m != 0){
            for(w=V; w>=0; w--){
                if (cur_N_ze[t][ww].A[w]> 0){
                    nnzp = cur_N_ze[t][ww].A[w];
                    break;
                }
            }
            if (nnzp > 0){
                break;
            }
        }
    }
    for (ww = 1; ww <= L; ww++){
        for (w = 1; w <= nz;w++){
            sum_1 += til[w][ww]*main_P_e[w][nzp];
        }
        sum_2[ww] = tjq0l_e[ww][nzp][nnzp] + sum_1;
        sum_1 = 0;
    }

    sum_1 = sum_2[1]; //поиск максимума
    for (ww=2; ww<=L; ww++){
        if (sum_2[ww] > sum_1){
            sum_1 = sum_2[ww];
        }
    }

    if(fabs(t_z - sum_1) <= 0.001 || (t_z - sum_1) > 0)
        return true;
    else return false;

}

//условие 19 для эффективной части алгоритма для z
bool check_condition_num(int t){ //номер 1 || 2
    //исследуется текущая группа z
    int w = 0, ww = 0, nzp = 0, nnzp = 0;
    int nz = 0;//число типов требований в текущей группе
    float sum_1 = 0, sum_2[L+1];
    //число партий в группе N_z
    for(ww=0; ww<Z; ww++){
        if (cur_N_z[t][ww].m != 0){
            nzp += (cur_N_z[t][ww].m); //решение строится с количеством партий m+1
            nz++;
        }
    }
    //количество требований входящих в последнюю партию
    for(ww=kz[t]; ww>=0; ww--){
        if (cur_N_z[t][ww].m != 0){
            for(w=V; w>=0; w--){
                if (cur_N_z[t][ww].A[w]> 0){
                    nnzp = cur_N_z[t][ww].A[w];
                    break;
                }
            }
            if (nnzp > 0){
                break;
            }
        }
    }
    for (ww = 1; ww <= L; ww++){
        for (w = 1; w <= nz;w++){
            sum_1 += til[w][ww]*main_P[w][nzp];
        }
        sum_2[ww] = tjq0l[ww][nzp][nnzp] + sum_1;
        sum_1 = 0;
    }

    sum_1 = sum_2[1]; //поиск максимума
    for (ww=2; ww<=L; ww++){
        if (sum_2[ww] > sum_1){
            sum_1 = sum_2[ww];
        }
    }

    if(fabs(t_z - sum_1) <= 0.001 || (t_z - sum_1) > 0)
        return true;
    else return false;

}

void PartyCount(){
    Zz = 0;
    int w = 0, ww = 0;
    for (w = 0; w < Z; w++){
        for (ww = 0; ww < Z; ww++){
            if (N_z[w][ww].m != 0){
                Zz++;
                break;
            }
        }
    }
}

int partyInGroupCount(int x){
    int ww = 0;
    int sum_1 = 0;
    for (ww = 0; ww < Z; ww++){
        if (N_z[x][ww].m != 0){
            sum_1++;
            break;
        }
    }
    return sum_1;
}

//------------------------------------------------------------

//------------------------эффективное решение--------------------------
void effective_middle(){
    ofstream out_file;
    out_file.open("C:\\output\\func2.txt", ios::out);
    out_file.close();
    int w = 0;
    /*for (z_e = 0; z_e <Z; z_e++){
        for (z = 0; z <Z; z++){
            N_z_effective[z_e][z].m = 0;
            N_z_effective[z_e][z].type = 0;
            for (w = 0; w <V; w++){
                N_z_effective[z_e][z].A[w] = 0;
            }
        }
    }*/
    //инициализация
    int tmp_z = 0;
    z_e = 0; //z' = 1;    
    for (z = 0; z< Z; z++){
        if (partyInGroupCount(z)){
            z_e = z;
            break;
        }
    }
    G2 = 0;
    G1 = 0;
    P[1] = P[2] = 0;//пустое множество

    int step = 1;
    finish_em = false;
    while(true){
        switch(step){
            case 1:
                //выбор партии для обмена
                z = 0;
                tmp_z = z_e+1;
                for (w = tmp_z; w <Z; w++){
                    if (partyInGroupCount(w)){
                        z = w;
                        break;
                    }
                }
                if (z == 0){
                    //больше нет партий для обмена
                    // пытаемся обменяться с Q
                    step = 40;
                } else {
                    //z = z_e + 1;                    
                    step = 2;
                    WriteGroupSwitchLog(z, z_e, 0, 1, 0, -1 ,"");
                }
                break;
            case 2:
                form_I_z(z, z_e);
                step = e_step_2();
                WriteFunc2ToFile();
                break;
            case 3:
                step = e_step_3();
                break;
            case 4:
                r_p[1] = main_R_e[Ip[1][1]][Jp[1][1]];
                step = 5;
                break;
            case 5:
                step = e_step_5(z_e, 1);
                break;
            case 6:
                step = e_step_6(z_e);
                break;
            case 7:
                step = e_step_7(z_e);
                break;
            case 8:
                step = e_step_8();
                break;
            case 9:
                step = e_step_5(z_e, 2);
                break;
            case 10:
                step = e_step_10();
                break;
            case 11:
                step = e_step_11();
                break;
            case 12:
                step = e_step_12();
                break;
            case 13:
                step = e_step_13(z, 1);
                break;
            case 14:
                step = e_step_14(z);
                break;
            case 15:
                step = e_step_15();
                break;
            case 16:
                step = e_step_16();
                break;
            case 17:
                step = e_step_13(z, 2);
                break;
            case 18:
                step = e_step_18();
                break;
            case 19:
                step = e_step_19(1);
                break;
            case 20:
                step = e_step_20(1);
                break;
            case 21:
                step = e_step_21_22(1);
                break;
            //case 22: break;
            case 23:
                step = e_step_23(1);
                WriteFunc2ToFile();
                break;
            case 24:
                step = e_step_24_25_26(1);
                break;
            //case 25: break;
            //case 26: break;
            case 27:
                step = e_step_27();
                break;
            case 28:
                step = e_step_19(2);
                break;
            case 29:
                step = e_step_20(2);
                break;
            case 30:
                step = e_step_21_22(2);
                break;
            //case 31: break;
            case 32:
                step = e_step_23(2);
                WriteFunc2ToFile();
                break;
            case 33:
                step = e_step_24_25_26(2);
                break;
            //case 34: break;
            //case 35: break;
            case 36:
                if (!P[1] && !P[2]){
                    step = 39;
                } else {
                    step = 37;
                }
                break;
            case 37:
                step = e_step_37();
                break;
            case 38:
                step = e_step_38();
                break;
            case 39:
                step = e_step_39();
                break;
            case 40:
                step = e_step_40();
                break;
            case 41:
                step = e_step_41_42_43(z_e);
                break;
            //case 42: break;
            //case 43: break;
            case 44:
                step = e_step_44();
                break;
            case 45:
                step = e_step_45();
                break;
            case 46:
                step = e_step_46_47();
                break;
            //case 47: break;
            case 48:
                step = e_step_48_49();
                break;
            //case 49: break;
            case 50:
                step = e_step_50(1);
                break;
            case 51:
                step = e_step_51();
                break;
            case 52:
                step = e_step_52();
                break;
            case 53:
                step = e_step_53(1);
                WriteFunc2ToFile();
                break;
            case 54:
                step = e_step_54_55(1);
                break;
           // case 55: break;
            case 56:
                step = e_step_56();
                break;
            case 57:
                step = e_step_50(2);
                break;
            case 58:
                step = e_step_58();
                break;
            case 59:
                step = e_step_53(2);
                WriteFunc2ToFile();
                break;
            case 60:
                step = e_step_54_55(2);
                break;
            //case 61: break;
            case 62:
                if (!P[1] && !P[2]){
                    step = 66;
                } else {
                    step = 63;
                }
                break;
            case 63:
                step = e_step_63();
                break;
            case 64:
                step = e_step_64_65();
                break;
            //case 65: break;
            case 66:
                step = e_step_66_67_68();
                break;
            //case 67: break;
            //case 68: break;
            case 69:
                step = e_step_69();
                break;
            case 70:
                step = e_step_70();
                break;
            case 71:
                step = e_step_71();
                break;
            case 72:
                step = e_step_72_73();
                break;
            //case 73: break;
        }        
        if (finish_em){
            writeEffGroupLog();
            break;
        }
    }
}
//---------------------------------------------------------------
void WriteGroupSwitchLog(int z, int z_e, int num, int step, int time, int q ,char* str){

    int w = 0, ww = 0;
    ofstream out_file;
    out_file.open("C://output//effective_log.txt", ios::app);
    out_file << "Step " << step << ":\n";
    out_file << "Switch between gr " << z_e << " and gr " << z << "\n";
    if (z_e>=0 && !num){
        out_file << "Group " << z_e << " was: "<< "\n";
        for (w = 0; w <Z; w++){
            if (N_z[z_e][w].m != 0){
                out_file << N_z[z_e][w].type << "\t" << N_z[z_e][w].m << "\t";
                for (ww = 1 ; ww < V; ww++){
                    out_file << N_z[z_e][ww].A[ww] << "\t";
                }
                out_file << "\n";
            }
        }
    }
    if (z>=0 && !num){
        out_file << "Group " << z << " was: "<< "\n";
        for (w = 0; w <Z; w++){
            if (N_z[z][w].m != 0){
                out_file << N_z[z][w].type << "\t" << N_z[z][w].m << "\t";
                for (ww = 1 ; ww < V; ww++){
                    out_file << N_z[z][ww].A[ww] << "\t";
                }
                out_file << "\n";
            }
        }
    }

    if (z_e>=0 && num){
        out_file << "Group " << z_e << " become: "<< "\n";
        for (w = 0; w <Z; w++){
            if (cur_N_ze[num][w].m != 0){
                out_file << cur_N_ze[num][w].type << "\t" << cur_N_ze[num][w].m << "\t";
                for (ww = 1 ; ww < V; ww++){
                    out_file << cur_N_ze[num][w].A[ww] << "\t";
                }
                out_file << "\n";
            }
        }
    }

    if (z>=0 && num){
        out_file << "Group " << z << " become: "<< "\n";
        for (w = 0; w <Z; w++){
            if (cur_N_z[num][w].m != 0){
                out_file << cur_N_z[num][w].type << "\t" << cur_N_z[num][w].m << "\t";
                for (ww = 1 ; ww < V; ww++){
                    out_file << cur_N_z[num][w].A[ww] << "\t";
                }
                out_file << "\n";
            }
        }
    }

    if (time){
        out_file << "Max time = " << time << "\t";
        if (time < t_z){
            out_file << " LESS then t_z" << "\n";
        } else {
            out_file << " MORE then t_z" << "\n";
        }
    }

    if (q >= 0){
        out_file << "Swith with Q" << "\n";
        out_file << Q[q].type << "\t" << Q[q].m << "\t";
        for (ww = 1 ; ww < V; ww++){
            out_file << Q[q].A[ww] << "\t";
        }
        out_file << "\n";
    }

    if (str){
        out_file << str << "\n";
    }

    out_file.close();
}


//nz - число типов требований в партиях группы
//nzp - число партий в последовательностях пи (строки в tjq0l)
void readFromFile(int param, int e, int nz, int nzp){ // 0 - P, 1 - R, 2 - tjq0l
                                    //e=0 для матриц main_P, e=1 для матриц main_P_e
    ifstream in;
    int w =0, ww = 0, www = 0;
    int ROW = 0, COL = 0;
    bool three = false;
    char* filename;
    int temp;
    int line[255];
    switch(param){
        case 0:
            filename = "C://output//P.txt";
            ROW = nz; COL = nzp; //COL = Z;
            for (w = 0; w<=N; w++){
                for (ww = 0; ww<=V; ww++){
                    if (e)
                        main_P_e[w][ww] = 0;
                    else
                        main_P[w][ww] = 0;
                }
            }
            break;
        case 1:
            filename = "C://output//R.txt";
            ROW = nz; COL = nzp;
            for (w = 0; w<=N; w++){
                for (ww = 0; ww<=V; ww++){
                    if (e)
                        main_R_e[w][ww] = 0;
                    else
                        main_R[w][ww] = 0;
                }
            }
            break;
        case 2:
            filename = "C://output//time//tjq0l.txt";
            ROW = nzp; COL = V;
            three = true;
            for(www = 1; www<=L; www++){ // по приборам
                for (w = 1; w <= nzp; w++){
                    for (ww = 1; ww <= QQ; ww++){
                        if (e)
                            tjq0l_e[www][w][ww] = 0;
                        else
                            tjq0l[www][w][ww] = 0;
                    }
                }
            }
            break;
    }
    int t = 0;
    temp = 0;
    char semicolon;
    in.open(filename, ios::in);
    if (three){ // матрица трехмерная
         if (in.is_open()){
             while ( in.good() )
             {
                for(www = 1; www<=L; www++){ // по приборам
                    for (w = 1; w <= nzp; w++){
                        for (ww = 1; ww <= QQ; ww++){ //!!!!!! hardcode
                            if (e){
                                in >> tjq0l_e[www][w][ww];
                            } else {                               
                                in >> tjq0l[www][w][ww];
                            }
                        }
                    }                    
                    in >> semicolon; //считает точку с запятой
                }
             }
             in >> semicolon; //считает точку с запятой или конец файла
         }

    } else {
        if (in.is_open()){
            temp = 1;
            while ( in.good() && temp <= N)
            {
                bool find = false;
                //for (w = 1; w <= nz; w++){
                w = 1;
                //temp = 1;
                while (w <= nz){
                    for (ww = 1; ww <= nzp; ww++){
                        in >> line[ww];                        
                        if (line[ww] > 0){
                            find = true;                            
                        }
                    }
                    temp++;
                    if (find){
                        find = false;
                        if (e){
                            if (param == 1)
                                for (ww = 1; ww <= nzp; ww++){
                                    main_R_e[w][ww] = line[ww];
                                }
                            else
                                for (ww = 1; ww <= nzp; ww++){
                                    main_P_e[w][ww] = line[ww];
                                }
                        } else {
                            if (param == 1)
                                for (ww = 1; ww <= nzp; ww++){
                                    main_R[w][ww] = line[ww];
                                }
                            else
                                for (ww = 1; ww <= nzp; ww++){
                                    main_P[w][ww] = line[ww];
                                }
                        }
                        w++;
                    }
                    if (temp > N){
                        break;
                    }
                }
                //in >> semicolon;
            }
            in.close();

        }else {
            cout<<"Error";
        }

    }
}

//надо будет разделять 2 решения с нижнего уровня
int find_j12(){// формироваться будет по Piz (j21 - 1 в степени)
    int w = 0, ww = 0, www = 0;
    int nzp = 0, position = 0;
    float sum_1[3];
    float max_elem = 0;
    int nj = 0;
    //число партий в группе N_z
    for(ww=0; ww<Z; ww++){
        if (N_z[z][ww].m != 0){
            nzp += (N_z[z][ww].m); //решение строится с количеством партий m+1
        }
    }
    float elem[V];
    for (w = 0; w <V ; w++){
        elem[w] = 0;
    }
    //правая часть неравенства
    for (w = 2; w <= nzp; w++){ //в формуле написано с 1
        elem[w] = 0;
        // вычисляем nj-1
        nj = 0;
        for (ww = 1; ww <= N; ww++){
            nj += main_R[ww][w-1];
        }
        for (ww = 1; ww <= L; ww++){
            elem[w] += tjq0l[ww][w][1];
            sum_1[0] = 0;
            for (www = 1; www <= kz[z]; www++){
                sum_1[0] +=til[www][ww]*main_P[www][w-1];
            }
            sum_1[0] +=tjq0l[ww][w-1][nj];
            elem[w] -= sum_1[0];
        }
    }

    //поиск максимума
    max_elem = elem[2];
    position = 2;
    for (w = 3; w <= nzp; w++){
        if(elem[w] > max_elem){
            max_elem = elem[w];
            position = w;
        }
    }

    //левая часть
    sum_1[1] = 0;
    for (ww = 1; ww <= L; ww++){
        sum_1[1] += tjq0l[ww][1][1];
    }
    if (sum_1[1] > max_elem){
        return 1;
    } else {
        return position;
    }

}

int find_j11(){// формироваться будет по Piz'
    int w = 0, ww = 0, www = 0;
    int nzp = 0, position = 0;
    float sum_1[3];
    float max_elem = 0;
    int nj = 0;
    //число партий в группе N_z
    for(ww=0; ww<Z; ww++){
        if (N_z[z_e][ww].m != 0){
            nzp += (N_z[z_e][ww].m); //решение строится с количеством партий m+1
        }
    }
    float elem[V];
    for (w = 0; w <V ; w++){
        elem[w] = 0;
    }
    //правая часть неравенства
    for (w = 2; w <= nzp; w++){ //в формуле написано с 1
        elem[w] = 0;
        // вычисляем nj-1
        nj = 0;
        for (ww = 1; ww <= N; ww++){
            nj += main_R_e[ww][w-1];
        }
        for (ww = 1; ww <= L; ww++){
            elem[w] += tjq0l_e[ww][w][1];
            sum_1[0] = 0;
            for (www = 1; www <= kz[z_e]; www++){
                sum_1[0] +=til[www][ww]*main_P_e[www][w-1];
            }
            sum_1[0] +=tjq0l_e[ww][w-1][nj];
            elem[w] -= sum_1[0];
        }
    }

    //поиск максимума
    max_elem = elem[2];
    position = 2;
    for (w = 3; w <= nzp; w++){
        if(elem[w] > max_elem){
            max_elem = elem[w];
            position = w;
        }
    }

    //левая часть
    sum_1[1] = 0;
    for (ww = 1; ww <= L; ww++){
        sum_1[1] += tjq0l_e[ww][1][1];
    }
    if (sum_1[1] > max_elem){
        return 1;
    } else {
        return position;
    }

}

void form_I_z(int one, int two){
    int w = 0, ww = 1, www = 1;
    for (w = 1; w <= N; w++){
        I_z[w] = I_ze[w] = -1;
    }
    for (w = 0; w < Z; w++){
        if(N_z[one][w].m != 0){
            I_z[ww] = N_z[one][w].type;
            ww++;
        }
        if(N_z[two][w].m != 0){
            I_ze[www] = N_z[two][w].type;
            www++;
        }
    }
}

int find_j21(int t){ //номер группы (j21 - 2 в степени)
    int w = 0, ww = 0, www = 0, wwww = 0;
    int nzp = 0, position = 0;
    int rj[V];
    for (w = 0; w< V; w++){
        rj[w] = 0;
    }
    float sum_1[3];
    float max_elem = 0;
    int nj = 0;
    //число партий в группе N_z
    for(ww=0; ww<Z; ww++){
        if (N_z[t][ww].m != 0){
            nzp += (N_z[t][ww].m); //решение строится с количеством партий m+1
        }
    }
    float elem[V];
    for (w = 0; w <V ; w++){
        elem[w] = 0;
    }
    for (w = 1; w <= nzp; w++){
        for(ww = 1; ww <= N; ww++){
            if (main_R_e[ww][w] > 0){
                rj[w] = main_R_e[ww][w] - 1;
            }
        }
    }

    for (w=1; w<=nzp; w++){ //по j
        elem[w] = 0;
        for (ww=1 ;ww<=L; ww++){ //по приборам
            for (www=1 ;www<=rj[w]; www++){
                sum_1[0] = 0;
                for (wwww=1 ;wwww<=kz[t]; wwww++){
                    sum_1[0] += til[wwww][ww]*main_P_e[wwww][w];
                }
                sum_1[1] = tjq0l_e[ww][w][www] - sum_1[0];
                sum_1[0] = tjq0l_e[ww][w][www-1] - sum_1[1];
            }
            elem[w] += sum_1[0];
        }
    }

    //максимум
    max_elem = elem[1];
    position = 1;
    for (w=2; w<=nzp; w++){
        if (elem[w] > max_elem){
            max_elem = elem[w];
            position = w;
        }
    }
    return position;
}

int find_j22(int t){ //номер группы
    int w = 0, ww = 0, www = 0, wwww = 0;
    int nzp = 0, position = 0;
    int rj[V];
    float sum_1[3];
    float max_elem = 0;
    int nj = 0;
    //число партий в группе N_z
    for(ww=0; ww<Z; ww++){
        if (N_z[t][ww].m != 0){
            nzp += (N_z[t][ww].m); //решение строится с количеством партий m+1
        }
    }
    float elem[V];
    for (w = 1; w <= nzp; w++){
        for(ww = 1; ww <= N; ww++){
            if (main_R[ww][w] != 0){
                rj[w] = main_R[ww][w] - 1;
            }
        }
    }

    for (w=1; w<=nzp; w++){ //по j
        elem[w] = 0;
        for (ww=1 ;ww<=L; ww++){ //по приборам
            for (www=1 ;www<=rj[w]; www++){
                sum_1[0] = 0;
                for (wwww=1 ;wwww<=kz[t]; wwww++){
                    sum_1[0] += til[wwww][ww]*main_P[wwww][w];
                }
                sum_1[1] = tjq0l[ww][w][www] - sum_1[0];
                sum_1[0] = tjq0l[ww][w][www-1] - sum_1[1];
            }
            elem[w] += sum_1[0];
        }
    }

    //максимум
    max_elem = elem[1];
    position = 1;
    for (w=2; w<=nzp; w++){
        if (elem[w] > max_elem){
            max_elem = elem[w];
            position = w;
        }
    }
    return position;
}

float find_f2(int nzp, int nnzp, int nj[Z], int nz){
    int w = 0, ww = 0, www = 0, wwww = 0;
    float sum_1[5];

    sum_1[0] = 0;
    //первая сумма
    for (w = 2; w <= L; w++){ //по приборам
        sum_1[0] += tjq0l[w][1][1];
    }

    //вторая сумма
    sum_1[1] = 0;
    for (w = 1; w <= L; w++){ //по приборам
        for (ww = 2; ww<= nzp; ww++){ //по j
            sum_1[2] = 0;
            for (www = 1; www <= nz; www++){
                sum_1[2] += til[www][w]*main_P[www][ww-1];
            }
           // sum_1[1] += sum_1[2]+tjq0l[w][ww-1][nj[ww-1]];
            sum_1[3] = tjq0l[w][ww][1] - sum_1[2] - tjq0l[w][ww-1][nj[ww-1]];
            sum_1[1] += sum_1[3];
        }
    }

    //третья сумма
    sum_1[2] = 0;
     for (w = 2; w <= L; w++){ //по приборам
         for (ww = 1; ww<= nzp; ww++){ //по j
             for (www = 2; www <= nj[ww]; www++){ //по q
                 sum_1[3] = 0;
                 for (wwww = 1; wwww <= nz; wwww++){ //по h
                     sum_1[3] +=til[wwww][w]*main_P[wwww][ww];
                 }
                 sum_1[3] = sum_1[3]+tjq0l[w][ww][www-1];
                 //sum_1[2] += sum_1[3]+tjq0l[w][ww][www-1];
                 sum_1[2] += tjq0l[w][ww][www] - sum_1[3];
             }
         }
     }

     //четвертая сумма
     sum_1[3] = 0;
     for (w = 1; w <= L; w++){ //по приборам
         sum_1[4] = 0;
         for (ww = 1; ww <= nz; ww++){ //по h
             sum_1[4] +=til[ww][w]*main_P[ww][nzp];
         }
         sum_1[4] +=tjq0l[w][nzp][nnzp];
         sum_1[3] += t_z - sum_1[4];
     }

     sum_1[4] = 0;
     for (w = 0; w<=3; w++){
         sum_1[4] += sum_1[w];
     }
     return sum_1[4];
}

float find_f2_e(int nzp, int nnzp, int nj[Z], int nz){
    int w = 0, ww = 0, www = 0, wwww = 0;
    float sum_1[5];

    sum_1[0] = 0;
    //первая сумма
    for (w = 2; w <= L; w++){ //по приборам
        sum_1[0] += tjq0l_e[w][1][1];
    }

    //вторая сумма
    sum_1[1] = 0;
    for (w = 1; w <= L; w++){ //по приборам
        for (ww = 2; ww<= nzp; ww++){ //по j
            sum_1[2] = 0;
            for (www = 1; www <= nz; www++){
                sum_1[2] += til[www][w]*main_P_e[www][ww-1];
            }
            sum_1[3] = tjq0l_e[w][ww][1] - sum_1[2] - tjq0l_e[w][ww-1][nj[ww-1]];
            sum_1[1] += sum_1[3];
            //sum_1[1] += sum_1[2]+tjq0l_e[w][ww-1][nj[ww-1]];
           // sum_1[1] = tjq0l_e[w][ww][1] - sum_1[1];
        }
    }

    //третья сумма
    sum_1[2] = 0;
     for (w = 2; w <= L; w++){ //по приборам
         for (ww = 1; ww<= nzp; ww++){ //по j
             for (www = 2; www <= nj[ww]; www++){ //по q
                 sum_1[3] = 0;
                 for (wwww = 1; wwww <= nz; wwww++){ //по h
                     sum_1[3] +=til[wwww][w]*main_P_e[wwww][ww];
                 }
                 sum_1[3] = sum_1[3]+tjq0l_e[w][ww][www-1];
                 //sum_1[2] += sum_1[3]+tjq0l[w][ww][www-1];
                 sum_1[2] += tjq0l_e[w][ww][www] - sum_1[3];
                 //sum_1[2] += sum_1[3]+tjq0l_e[w][ww][www-1];
                // sum_1[2] = tjq0l_e[w][ww][www] - sum_1[2];
             }
         }
     }

     //четвертая сумма
     sum_1[3] = 0;
     for (w = 1; w <= L; w++){ //по приборам
         sum_1[4] = 0;
         for (ww = 1; ww <= nz; ww++){ //по h
             sum_1[4] +=til[ww][w]*main_P_e[ww][nzp];
         }
         sum_1[4] +=tjq0l_e[w][nzp][nnzp];
         sum_1[3] += t_z - sum_1[4];
     }

     sum_1[4] = 0;
     for (w = 0; w<=3; w++){
         sum_1[4] += sum_1[w];
     }
     return sum_1[4];
}

int e_step_2(){
    int nz=0, nzp=0, nnzp=0, nj[Z];
    int w = 0, ww = 0;
    //число партий в группе N_z1
    for(ww=0; ww<Z; ww++){
        nj[ww] = 0;
        if (N_z[z][ww].m != 0){
            nzp += (N_z[z][ww].m); //решение строится с количеством партий m+1
        }
    }

    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (N_z[z][ww].m != 0){
            nz++;
        }
    }
    //количество требований входящих в последнюю партию
    for(ww=kz[z]; ww>=0; ww--){
        if (N_z[z][ww].m != 0){
            for(w=V; w>=0; w--){
                if (N_z[z][ww].A[w]> 0){
                    nnzp = N_z[z][ww].A[w];
                    break;
                }
            }
            if (nnzp > 0){
                break;
            }
        }
    }

    //передача на нижний уровень Nz и Nz'
    //вычисление критериев
    int col = 0;
    col = writeToFile(route_middle,z,0);
    low_level(col);
    //сначала передали на нижний уровень Nz
    //вычислили f2
    readFromFile(0,0,nz, nzp);
    readFromFile(1,0,nz, nzp);
    readFromFile(2,0,nz, nzp);

    //int nj[nzp+1];
    for (ww = 1; ww <= nzp+1; ww++){
        for (w = 1; w <= N; w++){
            if (main_R[w][ww] != 0){
                nj[ww] = main_R[w][ww];
                break;
            }
        }
    }

    f2_z[0] = find_f2(nzp,nnzp,nj,nz);

    //передали Nz'    
    col = writeToFile(route_middle,z_e,0);
    low_level(col);
    //получили соответствующие матрицы
    //вычисление критериев
    nz= 0, nzp = 0,nnzp=0;
    //число партий в группе N_z1
    for(ww=0; ww<Z; ww++){
        if (N_z[z_e][ww].m != 0){
            nzp += (N_z[z_e][ww].m); //решение строится с количеством партий m+1
        }
    }

    //количество требований входящих в последнюю партию
    for(ww=kz[z_e]; ww>=0; ww--){
        if (N_z[z_e][ww].m != 0){
            for(w=V; w>=0; w--){
                if (N_z[z_e][ww].A[w]> 0){
                    nnzp = N_z[z_e][ww].A[w];
                    break;
                }
            }
            if (nnzp > 0)
                break;
        }
    }

    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (N_z[z_e][ww].m != 0){
            nz++;
        }
    }
    readFromFile(0,1,nz, nzp);
    readFromFile(1,1,nz, nzp);
    readFromFile(2,1,nz, nzp);
    for (ww = 1; ww <= nzp+1; ww++){
        for (w = 1; w <= N; w++){
            if (main_R_e[w][ww] != 0){
                nj[ww] = main_R_e[w][ww];
                break;
            }
        }
    }
    Jp[1][1] = find_j11();
    f2_ze[0] = find_f2_e(nzp,nnzp,nj,nz);
    //вычислили f2'
    return 3;
}


int e_step_3(){
    //определение i11
    int w = 1;
    while (I_ze[w] != -1){
        if (main_P_e[w][Jp[1][1]] == 1){
            Ip[1][1] = w;
            i_p[1] = I_ze[w];
            break;
        }
        w++;
    }
    return 4;
}

//x = 1 (step 5), x = 2 (step 9)
int e_step_5(int t, int x){ //номер партии
    int w = 0, ww = 0;
    for (w = 0; w < Z; w++){
        if (N_z[t][w].m != 0 && N_z[t][w].type == i_p[x]){
            N_e[x].m = N_z[t][w].m;
            N_e[x].type = N_z[t][w].type;
            for (ww = 0; ww < V; ww++){
                N_e[x].A[ww] = N_z[t][w].A[ww];
            }
            break;
        }
    }
    if (x == 1)
        return 6;
    else
        return 10;
}

int e_step_6(int t){ //номер партии
    Jp[2][1] = find_j21(t);
    return 7;
}

int e_step_7(int t){ //номер партии
    //определение i21 (2 в степени)
    int w = 1;
    while (I_ze[w] != -1){
        if (main_P_e[w][Jp[2][1]] == 1){
            Ip[2][1] = w;
            i_p[2] = I_ze[w];
            break;
        }
        w++;
    }
    return 8;
}

int e_step_8(){
    //определение r'2
    int w = 0;
    r_p[2] = main_R_e[Ip[2][1]][Jp[2][1]];
    return 9;
}

int e_step_10(){ //номер партии
    Jp[1][2] = find_j12();
    return 11;
}

int e_step_11(){
    //определение i12 (1 в степени)
    int w = 1;
    while (I_z[w] != -1){
        if (main_P[w][Jp[1][2]] == 1){
            Ip[1][2] = w;
            ip[1] = I_z[w];
            break;
        }
        w++;
    }
    return 12;
}

int e_step_12(){
    //определение r1
    rp[1] = main_R[Ip[1][2]][Jp[1][2]];
    return 13;
}

//x = 1 (step 13), x = 2 (step 17)
int e_step_13(int t, int x){ //номер партии
    int w = 0, ww = 0;
    for (w = 0; w < Z; w++){
        if (N_z[t][w].m != 0 && N_z[t][w].type == ip[x]){
            N_num[x].m = N_z[t][w].m;
            N_num[x].type = N_z[t][w].type;
            for (ww = 0; ww < V; ww++){
                N_num[x].A[ww] = N_z[t][w].A[ww];
            }
            break;
        }
    }
    if (x == 1)
        return 14;
    else
        return 18;
}

int e_step_14(int t){
    Jp[2][2] = find_j22(t);
    return 15;
}

int e_step_15(){
    //определение i22
    int w = 1;
    while (I_z[w] != -1){
        if (main_P[w][Jp[2][2]] == 1){
            Ip[2][2] = w;
            ip[2] = I_z[w];
            break;
        }
        w++;
    }
    return 16;
}

int e_step_16(){
    //определение r1
    rp[2] = main_R[Ip[2][2]][Jp[2][2]];
    return 17;
}

int e_step_18(){
    if (ip[1] == i_p[1] && r_p[1] == rp[1]){
        //обмен не реализуется
        return 27;
    } else {
        return 19;
    }
}

int e_step_19(int num){ //num = 1 || num = 2(step=28)
    int w = 0, ww = 0;
    int position;
    k_zp[num] = kz[z_e];
    kzp[num] = kz[z];
    for (w = 0; w <Z ;w++){
        cur_N_ze[num][w].type = N_z[z_e][w].type;
        cur_N_ze[num][w].m = N_z[z_e][w].m;
        cur_N_z[num][w].type = N_z[z][w].type;
        cur_N_z[num][w].m = N_z[z][w].m;
        for (ww = 0; ww<V; ww++){
            cur_N_ze[num][w].A[ww] = N_z[z_e][w].A[ww];
            cur_N_z[num][w].A[ww] = N_z[z][w].A[ww];
        }
    }
    //для N1'
    for (w = 1; w <= N_e[num].m; w++){
        if (N_e[num].A[w] == r_p[num]){
            N_e[num].A[w] = 0;
            position = w;
            break;
        } //иначе ничего не изменяется
    }
    h_e = position;
    for (w = h_e; w <= N_e[num].m; w++){
        N_e[num].A[w] = N_e[num].A[w+1];
    }

    N_e[num].m--;
    if (N_e[num].m == 0){
        k_zp[num] --;
        //это партия исключается
        for (w = 0; w <Z ;w++){
            if(cur_N_ze[num][w].type == i_p[num] && cur_N_ze[num][w].m != 0){
                cur_N_ze[num][w].m = 0;
                break;
            }
        }
    } else {//иначе N1z' - это промежуточное решение (cur_N_ze[num])
        for (w = 0; w <Z ;w++){
            if(cur_N_ze[num][w].type == i_p[num] && cur_N_ze[num][w].m != 0){
                cur_N_ze[num][w].m = N_e[num].m;
                for (ww = 1 ; ww <V; ww++){
                    cur_N_ze[num][w].A[ww] = N_e[num].A[ww];
                }
                break;
            }
        }
    }

    //для N1
    for (w = 1; w <= N_num[num].m; w++){
        if (N_num[num].A[w] == rp[num]){
            N_num[num].A[w] = 0;
            position = w;
            break;
        } //иначе ничего не изменяется
    }
    h_e = position;
    for (w = h_e; w <= N_num[num].m; w++){
        N_num[num].A[w] = N_num[num].A[w+1];
    }

    N_num[num].m--;
    if (N_num[num].m == 0){
        kzp[num] --;
        //эта партия исключается
        for (w = 0; w <Z ;w++){
            if(cur_N_z[num][w].type == ip[num] && cur_N_z[num][w].m != 0){
                cur_N_z[num][w].m = 0;
                break;
            }
        }
    } else {//иначе N1z' - это промежуточное решение (cur_N_ze[num])
        for (w = 0; w <Z ;w++){
            if(cur_N_z[num][w].type == ip[num] && cur_N_z[num][w].m != 0){
                cur_N_z[num][w].m = N_num[num].m;
                for (ww = 1 ; ww <V; ww++){
                    cur_N_z[num][w].A[ww] = N_num[num].A[ww];
                }
                break;
            }
        }
    }

    //для i_p[1]
    position = 0;
    for (w = 0; w <= kzp[num]; w++){
        if (cur_N_z[num][w].type == i_p[num] && cur_N_z[num][w].m != 0){
            cur_N_z[num][w].m++;
            cur_N_z[num][w].A[cur_N_z[num][w].m] = r_p[num];
            position = 1;
        }
    }

    if (!position){ //такой тип не встречался
        //добавляем
        position = 0;
        for (w = 0; w <= kzp[num]; w++){
           if (cur_N_z[num][w].m == 0){
               cur_N_z[num][w].m = 1;
               cur_N_z[num][w].type = i_p[num];
               cur_N_z[num][w].A[cur_N_z[num][w].m] = r_p[num];
               kzp[num]++;
               position = 1;
               break;
           }
        }
        if (!position){ // в середину не вставилось
            //добавляю в конец
            cur_N_z[num][kzp[num]].m = 1;
            cur_N_z[num][kzp[num]].type = i_p[num];
            cur_N_z[num][kzp[num]].A[cur_N_z[num][kzp[num]].m] = r_p[num];
            kzp[num]++;
        }
    }

    //для ip[1]
    position = 0;
    for (w = 0; w <= k_zp[num]; w++){
        if (cur_N_ze[num][w].type == ip[num] && cur_N_ze[num][w].m != 0){
            cur_N_ze[num][w].m++;
            cur_N_ze[num][w].A[cur_N_ze[num][w].m] = rp[num];
            position = 1;
        }
    }

    if (!position){ //такой тип не встречался
        //добавляем
        position = 0;
        for (w = 0; w < k_zp[num]; w++){
           if (cur_N_ze[num][w].m == 0){
               cur_N_ze[num][w].m = 1;
               cur_N_ze[num][w].type = ip[num];
               cur_N_ze[num][w].A[cur_N_ze[num][w].m] = rp[num];
               k_zp[num]++;
               position = 1;
               break;
           }
        }
        if (!position){ // в середину не вставилось
            //добавляю в конец
            cur_N_ze[num][k_zp[num]].m = 1;
            cur_N_ze[num][k_zp[num]].type = ip[num];
            cur_N_ze[num][k_zp[num]].A[cur_N_ze[num][k_zp[num]].m] = rp[num];
            k_zp[num]++;
        }
    }

    if (num == 1){
        return 20;
    } else {
        return 29;
    }
}

int e_step_20(int num){ //num = 1 || num = 2(step=29)
    int nz = 0, nzp = 0;
    int w = 0, ww = 0;
    //на нижний уровень передается cur_N_z[num] и cur_N_ze[num]
    int col = 0;
    col = writeToFile(route_middle,num,1);
    low_level(col);

    //число партий в группе cur_N_z[1]
    for(ww=0; ww<Z; ww++){
        if (cur_N_z[num][ww].m != 0){
            nzp += (cur_N_z[num][ww].m); //решение строится с количеством партий m+1
        }
    }

    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (cur_N_z[num][ww].m != 0){
            nz++;
        }
    }
    readFromFile(0,0,nz, nzp);
    readFromFile(1,0,nz, nzp);
    readFromFile(2,0,nz, nzp);

    col = writeToFile(route_middle,num,2);
    low_level(col);

    nz = 0, nzp = 0;
    //число партий в группе cur_N_ze[1]
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[num][ww].m != 0){
            nzp += (cur_N_ze[num][ww].m); //решение строится с количеством партий m+1
        }
    }

    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[num][ww].m != 0){
            nz++;
        }
    }
    readFromFile(0,1,nz, nzp);
    readFromFile(1,1,nz, nzp);
    readFromFile(2,1,nz, nzp);

    if (num == 1){
        return 21;
    } else {
        return 30;
    }
}

int e_step_21_22(int num){ //num = 1 || num = 2(step=30-31)
    //проверка условий на длительность реализации для обоих расписаний
    //если не выполняется то переход к шагу 27 (36)
    bool y = false;
    y = check_condition_e(num);
    if (y){
        // проверка еще раз условия 19
        y = check_condition_num(num);
        if (y){
            if (num == 1)
                return 23;
            else
                return 32;
        } else {
            if (num == 1)
                return 27;
            else
                return 36;
        }

    } else{
        if (num == 1)
            return 27;
        else
            return 36;
    }
}

int e_step_23(int num){//num = 1 || num = 2(step=32)
    //расчет nz, nzp, nnzp
    int ww = 0, w = 0;
    int nzp = 0, nnzp = 0, nz = 0;
    int nj[Z];
    //число партий в группе N_z1
    for(ww=0; ww<Z; ww++){
        nj[ww] = 0;
        if (cur_N_z[num][ww].m != 0){
            nzp += (cur_N_z[num][ww].m); //решение строится с количеством партий m+1
        }
    }
    //количество требований входящих в последнюю партию
    for(ww=kz[z]; ww>=0; ww--){
        if (cur_N_z[num][ww].m != 0){
            for(w=V; w>=0; w--){
                if (cur_N_z[num][ww].A[w]> 0){
                    nnzp = cur_N_z[num][ww].A[w];
                    break;
                }
            }
            if (nnzp > 0)
                break;
        }
    }
    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (cur_N_z[num][ww].m != 0){
            nz++;
        }
    }
    //int nj[nzp+1];
    for (ww = 1; ww <= nzp+1; ww++){
        for (w = 1; w <= N; w++){
            if (main_R[w][ww] != 0){
                nj[ww] = main_R[w][ww];
                break;
            }
        }
    }

    fz2p[num] = find_f2(nzp,nnzp,nj,nz);


    nzp = 0, nnzp = 0, nz = 0;
    //число партий в группе N_z1'
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[num][ww].m != 0){
            nzp += (cur_N_ze[num][ww].m); //решение строится с количеством партий m+1
        }
    }
    //количество требований входящих в последнюю партию
    for(ww=kz[z]; ww>=0; ww--){
        if (cur_N_ze[num][ww].m != 0){
            for(w=V; w>=0; w--){
                if (cur_N_ze[num][ww].A[w]> 0){
                    nnzp = cur_N_ze[num][ww].A[w];
                    break;
                }
            }
            if (nnzp > 0)
                break;
        }
    }
    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[num][ww].m != 0){
            nz++;
        }
    }
    //int nj[nzp+1];
    for (ww = 1; ww <= nzp+1; ww++){
        for (w = 1; w <= N; w++){
            if (main_R_e[w][ww] != 0){
                nj[ww] = main_R_e[w][ww];
                break;
            }
        }
    }

    fz_2p[num] = find_f2_e(nzp,nnzp,nj,nz);
    //f_count++;

    if (num == 1)
        return 24;
    else
        return 33;
}

int e_step_24_25_26(int num){ //num = 1 || num = 2(step=33-34-35)
    int w = 0, ww = 0;
    leftGrad_ze[num] = fz_2p[num] - f2_ze[0];
    leftGrad_z[num] = fz2p[num] - f2_z[0];

    if (!(leftGrad_ze[num] < 0 && leftGrad_z[num] < 0)){
        P[num] = 0;//на всякий случай
        //return 27;
    } else {
        // решения рассматриваются как промежуточные
        //индекс p=1 включается в составл множества анализируемых решений
        P[num] = 1;
    }
    if (num == 1){
        return 27;
    } else {
        return 36;
    }
}

int e_step_27(){
    //два раза обмениваетмя тем же самым
    if(i_p[2] == i_p[1] && r_p[2] == r_p[1] && ip[2] == ip[1] && rp[2] == rp[1]){
        return 36;
    }
    if (i_p[2] == ip[2] && r_p[2] == rp[2]){           
        //обмен партиями не реализуется
        return 36;
    } else {
        return 28;
    }
}

int e_step_37(){
    //минимум градиента
    int p = 0, w = 0, ww = 0;
    if (leftGrad_ze[1] < leftGrad_ze[2]){
        p = 1;
    } else {
        p = 2;
    }
    //запоминаем локально решение
    for (w = 0; w < Z; w++){
        Nz_local[w].m = cur_N_z[p][w].m;
        Nz_e_local[w].m = cur_N_ze[p][w].m;
        Nz_local[w].type = cur_N_z[p][w].type;
        Nz_e_local[w].type = cur_N_ze[p][w].type;
        for (ww = 0; ww < V; ww++){
            Nz_e_local[w].A[ww] = cur_N_ze[p][w].A[ww];
            Nz_local[w].A[ww] = cur_N_z[p][w].A[ww];

        }
    }
    i_star = i_p[p];
    istar = ip[p];
    r_star = r_p[p];
    rstar = rp[p];

    return 38;
}

int e_step_38(){
    int min = 0;
    if (leftGrad_ze[1] < leftGrad_ze[2]){
        min = 1;
    } else {
        min = 2;
    }

    if (leftGrad_ze[min] < G1){
        //решение Nzp'* является более эффективным
        G1 = leftGrad_ze[min];
        l1 = z;
        i_ud[1] = i_star;
        i_dob[1] = istar;
        r_ud[1] = r_star;
        r_dob[1] = rstar;
    }
    return 39;
}

int e_step_39(){
    int w = 0; int ww = 0;
    bool find = false;
    ww = z+1;
    for (w = ww; w < Z; w++){
        //поиск следующей партии для обмена
        if (partyInGroupCount(w)){
            z = w;
            find = true;
            break;
        }
    }
    //z++;
    if (find){
        P[1] = P[2] = 0;
        return 2;
    } else {
        return 40;
    }
}

int e_step_40(){
    int w = 0, ww = 0;
    for (w = 0; w <Z; w++){
        cur_N_ze[1][w].m = N_z[z_e][w].m;
        cur_N_ze[1][w].type = N_z[z_e][w].type;
        cur_N_ze[2][w].m = N_z[z_e][w].m;
        cur_N_ze[2][w].type = N_z[z_e][w].type;
        for (ww = 0; ww <V; ww++){
            cur_N_ze[1][w].A[ww] = N_z[z_e][w].A[ww];
            cur_N_ze[2][w].A[ww] = N_z[z_e][w].A[ww];
        }
    }
    k_zp[1] = k_zp[2] = kz[z_e];
    P[1] = P[2] = 0;

    //передача Nz' на нижний уровень
    int nz=0, nzp=0, nnzp=0, nj[V];

    for (w = 0; w <V; w++){
        nj[w] = 0;
    }

    //число партий в группе N_z1
    for(ww=0; ww<Z; ww++){
        if (N_z[z_e][ww].m != 0){
            nzp += (N_z[z_e][ww].m); //решение строится с количеством партий m+1
        }
    }

    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (N_z[z_e][ww].m != 0){
            nz++;
        }
    }

    int col = 0;
    col = writeToFile(route_middle,z_e,0);
    low_level(col);
    //сначала передали на нижний уровень Nz
    //вычислили f2
    readFromFile(0,1,nz, nzp);
    readFromFile(1,1,nz, nzp);
    readFromFile(2,1,nz, nzp);

    form_I_z(0, z_e);

    //количество требований входящих в последнюю партию
    for(ww=kz[z_e]; ww>=0; ww--){
        if (N_z[z_e][ww].m != 0){
            for(w=V; w>=0; w--){
                if (N_z[z_e][ww].A[w]> 0){
                    nnzp = N_z[z_e][ww].A[w];
                    break;
                }
            }
            if (nnzp > 0)
                break;
        }
    }

    for (ww = 1; ww <= nzp+1; ww++){
        for (w = 1; w <= N; w++){
            if (main_R_e[w][ww] != 0){
                nj[ww] = main_R_e[w][ww];
                break;
            }
        }
    }
    f2_ze[0] = find_f2_e(nzp,nnzp,nj,nz);
    WriteFunc2ToFile();

    return 41;
}

int e_step_41_42_43(int t){ //номер партии
    Jp[1][1] = find_j11();
    Jp[2][1] = find_j21(t); //2 в степени

    //определение i21 (2 в степени)
    int w = 1;
    while (I_ze[w] != -1){
        if (main_P_e[w][Jp[2][1]] == 1){
            Ip[2][1] = w;
            i_p[2] = I_ze[w];
            break;
        }
        w++;
    }

    //определение i11 (1 в степени)
    while (I_ze[w] != -1){
        if (main_P_e[w][Jp[1][1]] == 1){
            Ip[1][1] = w;
            i_p[1] = I_ze[w];
            break;
        }
        w++;
    }

    r_p[1] = main_R_e[Ip[1][1]][Jp[1][1]];
    r_p[2] = main_R_e[Ip[2][1]][Jp[2][1]];
    return 44;
}

int e_step_44(){
    int tmp;
    //сформированы N'1 и N'2
    tmp = e_step_5(z_e,1);
    tmp = e_step_5(z_e,2);
    return 45;
}

int e_step_45(){
    //если Q пустое переход на шаг 72 (не указано в алгоритме)
    int w = 0;
    bool empty = true;

    for(w = 0; w< Z; w++){
        if (Q[w].m != 0){
            empty = false;
            break;
        }
    }
    if (empty){
        return 70;
    }
    kk = 0; //k = 1
    return 46;
}

int e_step_46_47(){
    hh = 1;
    ip[1] = Q[kk].type;
    return 48;
}

int e_step_48_49(){
    rp[1] = Q[kk].A[hh];
    if (i_p[1] == ip[1] && r_p[1] == rp[1]){
        //обмен не реализуется
        return 56;
    }
    return 50;
}

int e_step_50(int num){ //num = 1 || num = 2(step=57)
    int w = 0, ww = 0;
    int position;
    k_zp[num] = kz[z_e];
    kzp[num] = kz[z];
    for (w = 0; w <Z ;w++){
        cur_N_ze[num][w].type = N_z[z_e][w].type;
        cur_N_ze[num][w].m = N_z[z_e][w].m;
        //cur_N_z[num][w].type = N_z[z][w].type;
        //cur_N_z[num][w].m = N_z[z][w].m;
        for (ww = 0; ww<V; ww++){
            cur_N_ze[num][w].A[ww] = N_z[z_e][w].A[ww];
            //cur_N_z[num][w].A[ww] = N_z[z][w].A[ww];
        }
    }
    // поиск N1'
    for (w = 0 ; w < Z; w++){
        if(cur_N_ze[num][w].type == i_p[num] && cur_N_ze[num][w].m != 0){
            N_e[num].m = cur_N_ze[num][w].m;
            N_e[num].type = cur_N_ze[num][w].type;
            for (ww = 0 ; ww < V; ww++){
                N_e[num].A[ww] = cur_N_ze[num][w].A[ww];
            }
            break;
        }
    }
    //для N1'
    for (w = 1; w <= N_e[num].m; w++){
        if (N_e[num].A[w] == r_p[num]){
            N_e[num].A[w] = 0;
            position = w;
            break;
        } //иначе ничего не изменяется
    }
    h_e = position;
    for (w = h_e; w <= N_e[num].m; w++){
        N_e[num].A[w] = N_e[num].A[w+1];
    }

    N_e[num].m--;
    if (N_e[num].m == 0){
        k_zp[num] --;
        //это партия исключается
        for (w = 0; w <Z ;w++){
            if(cur_N_ze[num][w].type == i_p[num] && cur_N_ze[num][w].m != 0){
                cur_N_ze[num][w].m = 0;
                break;
            }
        }
    } else {//иначе N1z' - это промежуточное решение (cur_N_ze[num])
        for (w = 0; w <Z ;w++){
            if(cur_N_ze[num][w].type == i_p[num] && cur_N_ze[num][w].m != 0){
                cur_N_ze[num][w].m = N_e[num].m;
                for (ww = 1 ; ww <V; ww++){
                    cur_N_ze[num][w].A[ww] = N_e[num].A[ww];
                }
                break;
            }
        }
    }


    //для ip[1]
    position = 0;
    for (w = 0; w <= k_zp[num]; w++){
        if (cur_N_ze[num][w].type == ip[1] && cur_N_ze[num][w].m != 0){
            cur_N_ze[num][w].m++;
            cur_N_ze[num][w].A[cur_N_ze[num][w].m] = rp[1];
            position = 1;
        }
    }

    if (!position){ //такой тип не встречался
        //добавляем
        position = 0;
        for (w = 0; w <= k_zp[num]; w++){
           if (cur_N_ze[num][w].m == 0){
               cur_N_ze[num][w].m = 1;
               cur_N_ze[num][w].type = ip[1];
               cur_N_ze[num][w].A[cur_N_ze[num][w].m] = rp[1];
               k_zp[num]++;
               position = 1;
               break;
           }
        }
        if (!position){ // в середину не вставилось
            //добавляю в конец
            cur_N_ze[num][k_zp[num]].m = 1;
            cur_N_ze[num][k_zp[num]].type = ip[1];
            cur_N_ze[num][k_zp[num]].A[cur_N_ze[num][k_zp[num]].m] = rp[1];
            k_zp[num]++;
        }
    }

    if (num == 1){
        return 51;
    } else {
        return 58;
    }
}

int e_step_51(){    
    //решение Nz1'передается на 3 уровень для формирования расписания
    int col = 0;
    col = writeToFile(route_middle,1,2);
    low_level(col);
    int nz = 0, nzp = 0;
    int ww = 0;

    //число партий в группе cur_N_ze[1]
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[1][ww].m != 0){
            nzp += (cur_N_ze[1][ww].m); //решение строится с количеством партий m+1
        }
    }

    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[1][ww].m != 0){
            nz++;
        }
    }
    readFromFile(0,1,nz, nzp);
    readFromFile(1,1,nz, nzp);
    readFromFile(2,1,nz, nzp);
    return 52;
}

int e_step_52(){
    bool y = false;
    y = check_condition_e(1);
    if (y){
        return 53;
    } else {
        return 56;
    }
}

int e_step_53(int num){//num = 1 || num = 2(step=59)
    //расчет nz, nzp, nnzp
    int ww = 0, w = 0;
    int nzp = 0, nnzp = 0, nz = 0;

    //число партий в группе N_z1'
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[num][ww].m != 0){
            nzp += (cur_N_ze[num][ww].m); //решение строится с количеством партий m+1
        }
    }
    //количество требований входящих в последнюю партию
    for(ww=kz[z]; ww>=0; ww--){
        if (cur_N_ze[num][ww].m != 0){
            for(w=V; w>=0; w--){
                if (cur_N_ze[num][ww].A[w]> 0){
                    nnzp = cur_N_ze[num][ww].A[w];
                    break;
                }
            }
            if (nnzp > 0)
                break;
        }
    }
    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[num][ww].m != 0){
            nz++;
        }
    }
    int nj[V];
    for (ww = 1; ww <= nzp+1; ww++){
        for (w = 1; w <= N; w++){
            if (main_R_e[w][ww] != 0){
                nj[ww] = main_R_e[w][ww];
                break;
            }
        }
    }

    fz_2p[num] = find_f2_e(nzp,nnzp,nj,nz);
    //f_count++;

    if (num == 1)
        return 54;
    else
        return 60;
}

int e_step_54_55(int num){ //num = 1 || num = 2(step=60_61)
    int w = 0, ww = 0;
    leftGrad_ze[num] = fz_2p[num] - f2_ze[0];
    rightGrad_ze[num] = fz_2p[num] - f2_ze[0];
    if (leftGrad_ze[num] < 0){
        //полученное решение более эффективное чем исходное
        P[num] = num;
    }

    if (num == 1 && rightGrad_ze[num] >= 0){
        //решение не рассматривается в дальнейшем
    }
    if (num == 1)
        return 56;
    else
        return 62;

}

int e_step_56(){
    if (i_p[2] == ip[1] && r_p[2] == rp[1]){
       //обмен не выполняется
        return 62;
    } else {
        return 57;
    }
}

int e_step_58(){
    bool y = false;
    //перерача решение Nz2' на 3 уровень
    int col = 0;
    col = writeToFile(route_middle,2,2);
    low_level(col);
    int nz = 0, nzp = 0;
    int w=0, ww = 0;

    //число партий в группе cur_N_ze[2]
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[2][ww].m != 0){
            nzp += (cur_N_ze[2][ww].m); //решение строится с количеством партий m+1
        }
    }

    //число типов требований
    for(ww=0; ww<Z; ww++){
        if (cur_N_ze[2][ww].m != 0){
            nz++;
        }
    }
    readFromFile(0,1,nz, nzp);
    readFromFile(1,1,nz, nzp);
    readFromFile(2,1,nz, nzp);
    y = check_condition_e(2);
    if (y){
        return 59;
    } else {
        return 62;
    }
}

int e_step_63(){
    //минимум градиента
    int p = 0, w = 0, ww = 0;
    if (leftGrad_ze[1] < leftGrad_ze[2]){
        p = 1;
    } else {
        p = 2;
    }
    //запоминаем локально решение
    for (w = 0; w < Z; w++){
        //Nz_local[w].m = cur_N_z[p][w].m;
        Nz_e_local[w].m = cur_N_ze[p][w].m;
        //Nz_local[w].type = cur_N_z[p][w].type;
        Nz_e_local[w].type = cur_N_ze[p][w].type;
        for (ww = 0; ww < V; ww++){
            Nz_e_local[w].A[ww] = cur_N_ze[p][w].A[ww];
            //Nz_local[w].A[ww] = cur_N_z[p][w].A[ww];

        }
    }
    i_star = i_p[p];
    istar = ip[1];
    r_star = r_p[p];
    rstar = rp[1];

    return 64;
}

int e_step_64_65(){
    float min_grad = 0;
     if (leftGrad_ze[1] < leftGrad_ze[2]){
         min_grad = leftGrad_ze[1];
     } else {
         min_grad = leftGrad_ze[2];
     }
     if (min_grad < G2){
         G2 = min_grad;
         i_ud[2] = i_star;
         i_dob[2] = istar;
         r_ud[2] = r_star;
         r_dob[2] = rstar;
     } else {
         //текущее локально эффективное не является глобально эффективным
     }
     return 66;
}

int e_step_66_67_68(){
    hh++;
    P[1] = P[2] = 0;
    if (hh <= Q[kk].m){
        return 48;
    } else {
        kk++;
        if (kk<kq){
            return 46;
        } else {
            //получено локальное решение
        }
    }
    return 69;
}

int e_step_69(){
    int w = 0, ww = 0;
    //сравнение с точностью
    if (fabs(G1)<= epsilon && fabs(G2)<= epsilon){
    //получено эффективное решение по составу партий
        for (w = 0; w < Z; w++){
            N_z_effective[z_e][w].m = N_z[z_e][w].m;
            N_z_effective[z_e][w].type = N_z[z_e][w].type;
            for (ww = 1; ww < V; ww++){
                N_z_effective[z_e][w].A[ww] = N_z[z_e][w].A[ww];
            }
        }
        return 72;
    }
    return 70;
}

int e_step_70(){
    int w = 0, ww = 0;
    int position = 0, pos = 0;
    if (G1 == 0 && G2 == 0){
        return 72;
    }
    if (G1 < 0 || G2 < 0){
        if (G2 < G1){
            return 71;
        } else {
            for (w = 0; w<Z; w++){
                Nz_local[w].type = N_z[z][w].type;
                Nz_local[w].m = N_z[z][w].m;
                Nz_e_local[w].type = N_z[z_e][w].type;
                Nz_e_local[w].m = N_z[z_e][w].m;
                if (Nz_e_local[w].type == i_ud[1]){
                    pos = w;
                    N_e[1].type = Nz_e_local[w].type;
                    N_e[1].m = Nz_e_local[w].m;
                    for (ww=0; ww<V; ww++){
                        N_e[1].A[ww] = N_z[z_e][w].A[ww];
                    }
                }
                for (ww=0; ww<V; ww++){
                    Nz_local[w].A[ww] = N_z[z][w].A[ww];
                    Nz_e_local[w].A[ww] = N_z[z_e][w].A[ww];
                }
            }

            for (w = 0; w<=N_e[1].m; w++){
                if (N_e[1].A[w] == r_ud[1]){
                    N_e[1].A[w] = 0;
                    position = w;
                    break;
                }
            }
            for (w = position; w<=N_e[1].m; w++){
                N_e[1].A[w] = N_e[1].A[w+1];
            }
            N_e[1].m --;
            if (N_e[1].m != 0){
                z = l1;
                //сформировано промежуточное решение
            } else {
                Nz_e_local[pos].m = 0;
                kz[z_e]--;
            }
            //добавление

            position = 0;
            for (w = 0; w <= kz[z_e]; w++){
                if (Nz_e_local[w].type == i_dob[1] && Nz_e_local[w].m != 0){
                    Nz_e_local[w].m++;
                    Nz_e_local[w].A[Nz_e_local[w].m] = r_dob[1];
                    position = 1;
                    for (ww=0; ww<Z;ww++){
                       N_z[z_e][ww].type = Nz_e_local[ww].type;
                       N_z[z_e][ww].m = Nz_e_local[ww].m;
                       for (pos=0; pos <V; pos++){
                           N_z[z_e][ww].A[pos] = Nz_e_local[ww].A[pos];
                       }
                    }
                }
            }

            if (!position){ //такой тип не встречался
                //добавляем
                for (ww=0; ww<Z;ww++){
                   N_z[z_e][ww].type = Nz_e_local[ww].type;
                   N_z[z_e][ww].m = Nz_e_local[ww].m;
                   for (pos=0; pos <V; pos++){
                       N_z[z_e][ww].A[pos] = Nz_e_local[ww].A[pos];
                   }
                }
                position = 0;
                for (w = 0; w <= kz[z_e]; w++){
                   if ( N_z[z_e][w].m == 0){
                        N_z[z_e][w].m = 1;
                        N_z[z_e][w].type = i_dob[1];
                        N_z[z_e][w].A[N_z[z_e][w].m] = r_dob[1];
                       kz[z_e]++;
                       position = 1;
                       break;
                   }
                }
                if (!position){ // в середину не вставилось
                    //добавляю в конец
                    N_z[z_e][kz[z_e]].m = 1;
                    N_z[z_e][kz[z_e]].type = i_dob[1];
                    N_z[z_e][kz[z_e]].A[N_z[z_e][kz[z_e]].m] = r_dob[1];
                    kz[z_e]++;
                }
            }

            pos = 0;
            for (w = 0; w<Z; w++){
                if (Nz_local[w].type == i_dob[1] && Nz_local[w].m != 0){
                    pos = w;
                    N_num[1].type = Nz_local[w].type;
                    N_num[1].m = Nz_local[w].m;
                    for(ww=0; ww<V;ww++){
                        N_num[1].A[ww] = Nz_local[w].A[ww];
                    }
                    break;
                }
            }
            for (w = 0; w<=N_num[1].m; w++){
                if (N_num[1].A[w] == r_dob[1]){
                    N_num[1].A[w] = 0;
                    position = w;
                    break;
                }
            }
            for (w = position; w<=N_num[1].m; w++){
                N_num[1].A[w] = N_num[1].A[w+1];
            }
            N_num[1].m --;
            if (N_num[1].m != 0){
                //z = l1;
                //сформировано промежуточное решение
            } else {
                Nz_local[pos].m = 0;
                kz[z]--;
            }

            position = 0;
            for (w = 0; w <= kz[z_e]; w++){
                if (Nz_local[w].type == i_ud[1] && Nz_local[w].m != 0){
                    Nz_local[w].m++;
                    Nz_local[w].A[Nz_local[w].m] = r_ud[1];
                    position = 1;
                    for (ww=0; ww<Z;ww++){
                       N_z[z][ww].type = Nz_local[ww].type;
                       N_z[z][ww].m = Nz_local[ww].m;
                       for (pos=1; pos <V; pos++){
                           N_z[z][ww].A[pos] = Nz_local[ww].A[pos];
                       }
                    }
                }
            }

            if (!position){ //такой тип не встречался
                //добавляем
                for (ww=0; ww<Z;ww++){
                   N_z[z][ww].type = Nz_local[ww].type;
                   N_z[z][ww].m = Nz_local[ww].m;
                   for (pos=1; pos <V; pos++){
                       N_z[z][ww].A[pos] = Nz_local[ww].A[pos];
                   }
                }
                position = 0;
                for (w = 0; w <= kz[z]; w++){
                   if ( N_z[z][w].m == 0){
                        N_z[z][w].m = 1;
                        N_z[z][w].type = i_ud[1];
                        N_z[z][w].A[N_z[z][w].m] = r_ud[1];
                        kz[z]++;
                        position = 1;
                        break;
                   }
                }
                if (!position){ // в середину не вставилось
                    //добавляю в конец
                    N_z[z][kz[z]].m = 1;
                    N_z[z][kz[z]].type = i_ud[1];
                    N_z[z][kz[z]].A[N_z[z][kz[z]].m] = r_ud[1];
                    kz[z]++;
                }
            }
            //получено эффективное решение по составу партий
            for (w = 0; w < Z; w++){
                N_z_effective[z_e][w].m = N_z[z_e][w].m;
                N_z_effective[z][w].m = N_z[z][w].m;
                N_z_effective[z_e][w].type = N_z[z_e][w].type;
                N_z_effective[z][w].type = N_z[z][w].type;
                for (ww = 1; ww < V; ww++){
                    N_z_effective[z_e][w].A[ww] = N_z[z_e][w].A[ww];
                    N_z_effective[z][w].A[ww] = N_z[z][w].A[ww];
                }
            }
            G1 = 0;
            G2 = 0;
            P[1] = P[2] = 0;
            return 1;
        } //else
    }
}

int e_step_71(){
    int w = 0, ww = 0;
    int pos = 0, position = 0;
    for (w = 0; w < Z; w++){
       Nz_e_local[w].m = N_z[z_e][w].m;
       Nz_e_local[w].type = N_z[z_e][w].type;
       if (Nz_e_local[w].type == i_ud[2]){
           N_e[1].m = Nz_e_local[w].m;
           N_e[1].type = Nz_e_local[w].type;
           position = w;
           for (pos=0; pos<V;pos++){
               N_e[1].A[pos] =  Nz_e_local[w].A[pos];
           }
       }
       for (ww=0; ww<V;ww++){
            Nz_e_local[w].A[ww] = N_z[z_e][w].A[ww];
        }
    }

        for (w = 0; w<=N_e[1].m; w++){
            if (N_e[1].A[w] == r_ud[2]){
                N_e[1].A[w] = 0;
                pos = w;
                break;
            }
        }
        for (w = pos; w<=N_e[1].m; w++){
            N_e[1].A[w] = N_e[1].A[w+1];
        }
        N_e[1].m --;
        if (N_e[1].m != 0){
            //сформировано промежуточное решение
            Nz_e_local[position].m = N_e[1].m;
            for (w = 1; w <V ;w++){
                Nz_e_local[position].A[w] = N_e[1].A[w];
            }
        } else {
            Nz_e_local[position].m = 0;
            kz[z_e]--;
        }

        position = 0;
        for (w = 0; w <= kz[z_e]; w++){
            if (Nz_e_local[w].type == i_dob[2] && Nz_e_local[w].m != 0){
                Nz_e_local[w].m++;
                Nz_e_local[w].A[Nz_e_local[w].m] = r_dob[2];
                position = 1;
                for (ww=0; ww<Z;ww++){
                   N_z[z_e][ww].type = Nz_e_local[ww].type;
                   N_z[z_e][ww].m = Nz_e_local[ww].m;
                   for (pos=0; pos <V; pos++){
                       N_z[z_e][pos].m = Nz_e_local[pos].A[pos];
                   }
                }
            }
        }

        if (!position){ //такой тип не встречался
            //добавляем
            for (ww=0; ww<Z;ww++){
               N_z[z][ww].type = Nz_local[ww].type;
               N_z[z][ww].m = Nz_local[ww].m;
               for (pos=0; pos <V; pos++){
                   N_z[z][pos].m = Nz_local[pos].A[pos];
               }
            }
            position = 0;
            for (w = 0; w <= kz[z]; w++){
               if ( N_z[z][w].m == 0){
                    N_z[z][w].m = 1;
                    N_z[z][w].type = i_dob[2];
                    N_z[z][w].A[N_z[z][w].m] = r_dob[2];
                    kz[z]++;
                    position = 1;
                    break;
               }
            }
            if (!position){ // в середину не вставилось
                //добавляю в конец
                N_z[z][kz[z]].m = 1;
                N_z[z][kz[z]].type = i_dob[2];
                N_z[z][kz[z]].A[N_z[z][kz[z]].m] = r_dob[2];
                kz[z_e]++;
            }
        }

        for (w=0;w<Z;w++){
            if (Q[w].m != 0 && Q[w].type == i_dob[2]){
                pos = w; //запомнили элемент
                break;
            }
        }
        for (w = 0; w<=Q[pos].m; w++){
            if (Q[pos].A[w] == r_dob[2]){
                Q[pos].A[w] = 0;
                ww = w;
                break;
            }
        }
        for (w = ww; w<=Q[pos].m; w++){
            Q[pos].A[w] = Q[pos].A[w+1];
        }
        Q[pos].m --;
        if (Q[pos].m != 0){
            //сформировано промежуточное решение
        } else {
            Q[pos].m = 0;
            kq--;
        }

        position = 0;
        for (w = 0; w < kq; w++){
            if (Q[w].type == i_ud[2] && Q[w].m != 0){
                Q[w].m++;
                Q[w].A[Q[w].m] = r_ud[2];
                position = 1;
            }
        }

        if (!position){ //такой тип не встречался
            //добавляем
            position = 0;
            for (w = 0; w < kq; w++){
               if ( Q[w].m == 0){
                    Q[w].m = 1;
                    Q[w].type = i_ud[2];
                    Q[w].A[Q[w].m] = r_ud[2];
                    kq++;
                    position = 1;
                    break;
               }
            }
            if (!position){ // в середину не вставилось
                //добавляю в конец
                Q[kq].m = 1;
                Q[kq].type = i_ud[2];
                Q[kq].A[Q[kq].m] = r_ud[2];
                kq++;
            }
        }
        //получено эффективное решение по составу партий
        for (w = 0; w < Z; w++){
            N_z_effective[z_e][w].m = N_z[z_e][w].m;
            N_z_effective[z_e][w].type = N_z[z_e][w].type;
            for (ww = 1; ww < V; ww++){
                N_z_effective[z_e][w].A[ww] = N_z[z_e][w].A[ww];
            }
        }

        G1 = 0;
        G2 = 0;
        P[1] = P[2] = 0;
        return 1;
}

int e_step_72_73(){
    int w = 0;
    bool find = false;
    z_e++;
    for (w = z_e; w <Z; w++){
        if (partyInGroupCount(w)){
            z_e = w;
            find = true;
            break;
        }
    }
    //z_e++;
    if (find){
        return 1;
    } else {
        //получено эффективное решение
        finish_em = true;
    }
}

//-------------------------Нижний уровень------------------
void start_time (int **P, int **R, int n_z, int n_p_z)
/*Вычисление матрицы времени (t_jq_0l) начала обработки q-ых требований
партий, занимающих в последовательности Pi_l j-ю позицию, q=1..r_ij*/
{
    int q,i1,j1,n[n_p_z+1];
    int i,j,l,h;
    float sum1, sum2, sum3, sum4;
    float t_jq_0l[n_p_z+1][QQ+1][L+1];   //нач. обработки q-ых треб. партии, заним в Pi_l
                                        //j-ю позицию


    for (i=1; i<=QQ; i++)
    {
        for (j=1; j<=n_p_z; j++)
        {
            for (l=1; l<=L; l++)
            {
                t_jq_0l[j][i][l]=0;
            }
        }
    }
    l=1;

    //для первого прибора
    for (j=1;j<=n_p_z;j++)
    {
        n[j]=0;

        for (h=1; h<=n_z; h++)
        {
            n[j]+=R[h][j];
        }


        for (q=1; q<=n[j]; q++)
         {
                sum1=0;
                for (h=1;h<=n_z;h++)
                {
                    sum1+=t_ii_l[h][1]*P[h][1];
                }
                sum2=0;
                for (int f=1;f<=j-1;f++)
                for (h=1;h<=n_z;h++)
                {
                    sum2+=t_i_l[h][1]*R[h][f];
                }
                sum3=0;
                for (h=1;h<=j-1;h++)
                {
                    int ir1=0,ir2=0;
                    for (i=1;i<=n_z;i++)
                    {
                        if (P[i][h]==1)
                            ir1=i;
                        if (P[i][h+1]==1)
                            ir2=i;
                    }
                    if ((ir1!=0)&&(ir2!=0))
                    {
                        sum3+=t_ik_l[ir1][ir2][1];
                    }
                }
                sum4=0;
                for (h=1;h<=n_z;h++)
                {
                    sum4+=t_i_l[h][1]*P[h][j];
                }
                t_jq_0l[j][q][1]=sum1+sum2+sum3+(q-1)*sum4;
            }
        }

    //для остальных приборов
    for (l=2; l<=L; l++)
    {
        for (j=1;j<=n_p_z;j++)
        {
             for (i=1; i<=n_z; i++)
            {
                 if (P[i][j]==1)
                {
                     i1=i;
                }
            }
            n[j]=0;

            for (h=1; h<=n_z; h++)
            {
                n[j]+=R[h][j];
            }


            for (q=1; q<=n[j]; q++)
             {
                 //для первого требования в первой позиции
                 if ((q==1)&&(j==1))
                 {
                     sum1=0;
                     sum2=0;
                     for (h=1; h<=n_z; h++)
                     {
                         sum1+=t_i_l[h][l]*P[h][1];
                         sum2+=t_i_l[h][l-1]*P[h][1];
                     }
                     sum2+=t_jq_0l[1][1][l-1];
                     t_jq_0l[1][1][l]=std::max(sum1,sum2);
                 }
                 else
                 {
                     //для первого требования в непервой позиции
                     if ((q==1)&&(j!=1))
                    {
                        sum1=0;
                        sum2=0;
                        for (h=1; h<=n_z; h++)
                        {
                             sum1+=t_i_l[h][l]*P[h][j-1];
                             sum2+=t_i_l[h][l-1]*P[h][j];
                        }
                        int ir1=0,ir2=0;
                        for (i=1;i<=n_z;i++)
                        {
                            if (P[i][j-1]==1)
                                ir1=i;
                            if (P[i][j]==1)
                                ir2=i;
                        }
                        sum1+=t_jq_0l[j-1][n[j-1]][l];
                        sum1+=t_ik_l[ir1][ir2][l];
                        sum2+=t_jq_0l[j][1][l-1];
                        t_jq_0l[j][1][l]=max(sum1,sum2);
                    }
                     //все остальные требования
                     else
                        {
                         sum1=0;
                         sum2=0;
                         for (h=1; h<=n_z; h++)
                         {
                             sum1+=t_i_l[h][l]*P[h][j];
                             sum2+=t_i_l[h][l-1]*P[h][j];
                         }
                         sum1+=t_jq_0l[j][q-1][l];
                         sum2+=t_jq_0l[j][q][l-1];
                         t_jq_0l[j][q][l]=max(sum1,sum2);
                        }
                    }
             }
         }

    }

    //записать результаты в файл
    std::ofstream tjq0l (route_tjq0l);
    for (l=1; l<=L; l++)
    {
        for (j=1; j<=n_p_z; j++)
        {
            for (q=1; q<=QQ; q++)
            {
                if (t_jq_0l[j][q][l]>0)
                {
                    tjq0l << t_jq_0l[j][q][l] << "\t";
                   // std::cout << t_jq_0l[j][q][l] << " ";
                }
                else
                {
                    tjq0l << "0\t";
                    //std::cout << "0 ";
                }
            }
            tjq0l << "\n";
            //std::cout << std::endl;
        }
        tjq0l << ";\n";
        //std::cout << ";" << std::endl;
   }
    tjq0l.close();
}

float objective(int **P, int **R, int n_z, int n_p_z)
/*Вычисление целевой функции*/
{
    float sum1, sum2, sum3;
    int l, j, q, h, n[n_p_z+1];
    float t_jq_0l[n_p_z+1][QQ+1][L+1];

    std::ifstream tjq0l(route_tjq0l);
    char c;
    for (l=1; l<=L; l++)
    {
        for (j=1; j<=n_p_z; j++)
        {
            for (q=1; q<=QQ; q++)
            {
                tjq0l >> t_jq_0l[j][q][l];
            }
        }
        tjq0l >> c;
    }

    sum1=0;
    sum2=0;
    sum3=0;

    for (j=1; j<=n_p_z; j++)
    {
        for (h=1; h<=n_z; h++)
        {
            n[j]=R[h][j];
        }
    }

    for (l=1; l<=L; l++)
    {
        sum1+=t_jq_0l[1][1][l];
        for (j=2; j<=n_p_z; j++)
        {
            for (h=1; h<=n_z; h++)
            {
                sum2+=t_i_l[h][l]*P[h][j-1];
            }
            if (n[j-1]!=0)
            {
                sum2+=t_jq_0l[j-1][n[j-1]][l];
            }
            sum2=t_jq_0l[j][1][l]-sum2;
        }
    }

    for (l=2; l<=L; l++)
    {
        for (j=1; j<=n_p_z; j++)
        {
            for (q=2; q<=n[j]; q++)
            {
                for (h=1; h<=n_z; h++)
                {
                    sum3+=t_i_l[h][l]*P[h][j];
                }
                sum3+=t_jq_0l[j][q-1][l];
                sum3=t_jq_0l[j][q][l]-sum3;
            }
        }
    }
    return (sum1+sum2+sum3);
}

void start_part(int **P, int **R, int n_z, int n_p_z)
/*Вычисление матрицы времени (t_ji_nl)начала обработки i-тых партий,
 *находящихся в последовательности Pi_l на j-ой позиции*/
{
    float t_ji_nl [n_p_z+1][n_z+1][L+1],
            t_jq_0l[n_p_z+1][QQ+1][L+1];
    int l, j, q, i, h, f;
    float sum1, sum2, sum3;

    std::ifstream tjq0l(route_tjq0l);
    char c;
    for (l=1; l<=L; l++)
    {
        for (j=1; j<=n_p_z; j++)
        {
            for (q=1; q<=QQ; q++)
            {
                tjq0l >> t_jq_0l[j][q][l];
            }
        }
        tjq0l >> c;
    }

    for (i=1; i<=n_z; i++)
    {
        for (j=1; j<=n_p_z; j++)
        {
            for (l=1; l<=L; l++)
            {
                t_ji_nl[j][i][l]=0;
            }
        }
    }

    l=1;

    //для первого прибора
    for (j=1; j<=n_p_z; j++)
    {
        for (i=1; i<=n_z; i++)
        {
            sum1=0;
            sum2=0;
            sum3=0;

            for (h=1; h<=n_z; h++)
            {
                sum1+=t_ii_l[h][1]*P[h][1];
            }

            for (f=1; f<=j-1; f++)
            {
                for (h=1;h<=n_z; h++)
                {
                    sum2+=t_i_l[h][1]*R[h][f];
                }
            }

            for (h=1; h<=j-1; h++)
            {
                int ir1=0,ir2=0;
                for (i=1;i<=n_z;i++)
                {
                    if (P[i][h]==1)
                        ir1=i;
                    if (P[i][h+1]==1)
                        ir2=i;
                }
                if ((ir1!=0)&&(ir2!=0))
                {
                    sum3+=t_ik_l[ir1][ir2][1];
                }
            }

            t_ji_nl[j][i][1]=sum1+sum2+sum3;
        }
    }

    //для остальных приборов
    for (l=2; l<=L; l++)
    {
        for (j=1; j<=n_p_z; j++)
        {
            for (i=1; i<=n_z; i++)
            {
                if (j==1)
                {
                    //для партии в первой позиции
                    t_ji_nl[1][i][l]=t_jq_0l[1][1][l];
                }
                else
                {
                    //для всех остальных партий
                    t_ji_nl[j][i][l]=t_jq_0l[j][1][l];
                }
            }
        }
    }

    //записать результаты в файл
/*    std::ofstream tjinl (route_tjinl);
    for (l=1; l<=L; l++)
    {
        for (j=1; j<=n_p_z; j++)
        {
            for (i=1; i<=n_z; i++)
            {
                if (t_ji_nl[j][i][l]>0)
                {
                    tjinl << t_ji_nl[j][i][l] << "\t";
                    std::cout << t_ji_nl[j][i][l] << " ";
                }
                else
                {
                    tjinl << "0\t";
                    std::cout << "0 ";
                }
            }
            tjinl << "\n";
            std::cout << std::endl;
        }
        tjinl << ";\n";
        std::cout << ";" << std::endl;
   }
    tjinl.close();*/
}



//*****************************************************************
int low_level(int x)
{
    //QCoreApplication a(argc, argv);
    //--------------------Входные данные-------------------------
    //------------------N_z={[i,m_i,(A)_i]_k|k=1..k_z}-----------
    struct Group
    {
        int i,          //тип требования
            mi,         //количество партий
            *Ai;        //вектор количества требований в каждой партии
    };

    int k_z;    //число типов требований в группе N_z
    int *I_z;   //вектор типов требований, партии которых принадлежат группе N_z
    Group *N_z; //группы



    //-----------------------Переменные--------------------------
    int s,      //шаг алгоритма
        g,      //индекс текущего промежуточного решения
        n_p_z2, //количество партий требований, размещённых в последовательностях
                //Pi_l на предыдущих s шагах реализации алгоритма
        i2,     //индекс в массиве I_z
         v,     //индекс текущего столбца матриц P_z(S) и R_z(S)
        k,      //k=1..k_z
        v_max,
            h,  //индекс в векторе А [1..mi]
        ***P,
        ***R,
        **P_eff,
        **R_eff,
        **P_temp,
        **R_temp,
        n_z,
        n_p_z,
        step,
           l,
           mi,
        w_max=3,
        g_max;
    float *f3,
          f3_eff,
          *w;
    char iter[11],
            logfile[40];

    //логирование
    sprintf(iter,"%d", iteration);
    strcpy(logfile,route_log);
    strcat(strcpy(logfile,route_log),iter);
    strcat(logfile,".txt");
    ofstream log (logfile);
    //----------------Тестовые входные данные----------------------
    //файл должен быть записан правильно!!!!
    k_z=0;
    ifstream in (route_middle);
    /*char *str = new char [1024];
    while (!in.eof())
    {
        in.getline(str,1024,'\n');
        k_z++;
    }
    in.close();
    in.open (route_middle);
    */
    k_z = x;
    I_z=new int[k_z+1];
    N_z=new Group[k_z+1];
    k=0;
    while (k<k_z)
    {
        k++;
        in >> I_z[k];
        N_z[k].i=I_z[k];
        log << "i:" << I_z[k];
        in >> N_z[k].mi;
        log << " m:" << N_z[k].mi;
        N_z[k].Ai=new int[N_z[k].mi];
        for (int i=1; i<=N_z[k].mi; i++)
        {
            in >> N_z[k].Ai[i];
        }
        log << " A:";
        for (int i=1; i<=N_z[k].mi;i++)
        {
            log << " " << N_z[k].Ai[i];
        }
        log << std::endl;
    }
    in.close();
    log << "Number of types in group: " << k << std::endl;

    //---------------Входные данные (матрицы времён)---------------
    ifstream file_til (route_til);
    for (int i=1; i<=II; i++)
    {
        for (int l=1; l<=L; l++)
        {
            file_til >> t_i_l[i][l];
        }
    }
    file_til.close();

    char c;
    ifstream tikl (route_tikl);
    for (int l=1; l<=L; l++)
    {
        for (int i=1; i<=II; i++)
        {
            for (int k=1; k<=II; k++)
            {
                tikl >> t_ik_l[i][k][l];
            }
        }
        tikl >> c;
    }
    tikl.close();

    ifstream tiil (route_tiil);
    for (int i=1; i<=II; i++)
    {
        for (int l=1; l<=L; l++)
        {
            tiil >> t_ii_l[i][l];
        }
    }
    tiil.close();


    //------------Начальная инициализация параметров-----------------
    s=0;
    g=0;
    g_max=20;
    k=1;
    v_max=0;
    n_z=II;

    n_p_z=0;
    n_p_z2=0;
    for (int i=1;i<=k_z;i++)
    {
        n_p_z+=N_z[i].mi;
    }

    //выделение памяти для матриц P и R
    P=new int **[n_z+1];
    for (int i=0;i<n_z+1;i++)
    {
        P[i]=new int *[n_p_z+1];
        for (int j=0; j<n_p_z+1; j++)
        {
            P[i][j]=new int [g_max+1];
        }
    }

    R=new int **[n_z+1];
    for (int i=0;i<n_z+1;i++)
    {
        R[i]=new int *[n_p_z+1];
        for (int j=0; j<n_p_z+1; j++)
        {
            R[i][j]=new int [g_max+1];
        }
    }

    for (int i=1;i<=n_z;i++)
    {
        for (int j=1;j<=n_p_z;j++)
        {
            for (int g=0; g<=g_max; g++)
            {
                P[i][j][g]=R[i][j][g]=0;
            }

        }
    }

    P_eff=new int *[n_z+1];
    for (int i=0;i<n_z+1;i++)
    {
        P_eff[i]=new int[n_p_z+1];
    }

    R_eff=new int *[n_z+1];
    for (int i=0;i<n_z+1;i++)
    {
        R_eff[i]=new int[n_p_z+1];
    }

    for (int i=1;i<=n_z;i++)
        for (int j=1;j<=n_p_z;j++)
        {
            P_eff[i][j]=R_eff[i][j]=0;
        }

    log << std::endl << "P: " << std::endl;
    for (int i=1; i<=n_z; i++)
    {
        for (int j=1; j<=n_p_z; j++)
        {
            log << P_eff[i][j] << " ";
        }
        log << std::endl;
    }
    log << std::endl << "R: " << std::endl;
    for (int i=1; i<=n_z; i++)
    {
        for (int j=1; j<=n_p_z; j++)
        {
            log << R_eff[i][j] << " ";
        }
        log << std::endl;
    }

    P_temp=new int *[n_z+1];
    for (int i=0;i<n_z+1;i++)
    {
        P_temp[i]=new int[n_p_z+1];
    }

    R_temp=new int *[n_z+1];
    for (int i=0;i<n_z+1;i++)
    {
        R_temp[i]=new int[n_p_z+1];
    }

    for (int i=1;i<=n_z;i++)
    {
        for (int j=1;j<=n_p_z;j++)
        {
            P_temp[i][j]=R_temp[i][j]=0;
        }
    }

    f3=new float [g_max+1];
    w=new float [g_max+1];

    //--------------------------Алгоритм-------------------------

    k=1;
    int stop=0;
    step=1;
    //первый тип полностью заносится в расписание, так как от перестановки
    //партий требований одного типа целевая функция не изменяется
    for (int i1=1; i1<=N_z[k].mi; i1++)
    {
        v_max++;
        P_eff[N_z[k].i][v_max]=P[N_z[k].i][v_max][0]=1;
        R_eff[N_z[k].i][v_max]=R[N_z[k].i][v_max][0]=N_z[k].Ai[i1];
        start_time(P_eff, R_eff, n_z, n_p_z);
        f3_eff=objective(P_eff, R_eff, n_z, n_p_z);
        start_part(P_eff, R_eff, n_z, n_p_z);
        s++;
        n_p_z2++;
    }
    log << "\n Pervii tip trebovanii polnostiyu dobavlen: \n";
    log << "s: " << s << " g:" << g << "\n";
    log << std::endl << "P: " << std::endl;
    for (int i=1; i<=n_z; i++)
    {
        for (int j=1; j<=n_p_z; j++)
        {
            log << P_eff[i][j] << " ";
        }
        log << std::endl;
    }
    log << std::endl << "R: " << std::endl;
    for (int i=1; i<=n_z; i++)
    {
        for (int j=1; j<=n_p_z; j++)
        {
            log << R_eff[i][j] << " ";
        }
        log << std::endl;
    }
    log << "Kriterii (f3): " << f3_eff << "\n";

    k=2;
    while ((stop!=1)&&(k<=k_z))
    {
        switch (step)
        {
            case 1:
            {
                v_max++;
                mi=N_z[k].mi;
                h=1;
                step=2;
                break;
            }
            case 2:
            {
                g=0;
                for (int i=1;i<=n_z; i++)
                {
                    if (i==I_z[k])
                    {
                        P[i][v_max][g]=1;
                        R[i][v_max][g]=N_z[k].Ai[h];
                    }
                    else
                    {
                        P[i][v_max][g]=0;
                        R[i][v_max][g]=0;
                    }

                }
                step=3;

                //проверка
                log << "Dobavlena partiya \n";
                log << "P:" << std::endl;
                for (int i=1;i<=n_z;i++)
                {
                    for (int j=1; j<=n_p_z;j++)
                    {
                        log << P[i][j][g] << " ";
                    }
                    log << std::endl;
                }
                log << std::endl;
                log << "R:" << std::endl;
                for (int i=1;i<=n_z;i++)
                {
                    for (int j=1; j<=n_p_z;j++)
                    {
                        log << R[i][j][g]<< " ";
                    }
                    log << std::endl;
                }
                break;
            }
            case 3:
            {
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        P_temp[i][j]=P[i][j][g];
                        R_temp[i][j]=R[i][j][g];
                    }
                }
                //вычисление элементов матрицы t_jq_0l c l=1..L
                start_time(P_temp,R_temp,n_z,n_p_z);
                step=4;
                break;
            }
            case 4:
            {
                //определение критерия f(s)
                f3[g]=objective(P_temp,R_temp,n_z,n_p_z);
                log << std::endl;
                log << "Kriterii (f3): " << f3[g] << std::endl;

                //вычисление элементов матрицы t_ij_nl
                start_part(P_temp,R_temp, n_z, n_p_z);

                //фиксация локально эффективного решения
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        P_eff[i][j]=P[i][j][g];
                        R_eff[i][j]=R[i][j][g];
                    }
                }
                f3_eff=f3[g];
                //фиксируется локально-эффективное решение
                log << "--------Lokalno-effektivnoe reshenie: s=" << s;
                ofstream P_file(route_P);
                log << " g=" << g << std::endl << "P: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        log << P_eff[i][j] << " ";
                        P_file << P_eff[i][j] << "\t";
                    }
                    log << std::endl;
                    P_file << "\n";
                }
                log << std::endl;
                ofstream R_file (route_R);
                log << "R: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        log << R_eff[i][j] << " ";
                        R_file << R_eff[i][j] << "\t";
                    }
                    log << std::endl;
                    R_file << "\n";
                }
                log << std::endl;
                step=5;
                break;
            }
            case 5:
            {
                if (v_max==1)
                {
                    //если добавляемая в поледовательность Pi_l партия 1-ая
                    step=13;
                }
                else
                {
                    step=6;
                }
                break;
            }
            case 6:
            {
                g=1;
                v=v_max;
                step=7;
                break;
            }
            case 7:
            {
                //изменение порядка партий
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        P[i][j][g]=P[i][j][g-1];
                        R[i][j][g]=R[i][j][g-1];
                    }
                }

                int i1;
                for (int i=1; i<=n_z; i++)
                {
                    if (P[i][v-1][g-1]==1)
                    {
                        i1=i;
                    }
                }

                //std::cout << N_z[k].i << std::endl;
                P[i1][v-1][g]=0;
                P[N_z[k].i][v][g]=0;
                P[i1][v][g]=1;
                R[N_z[k].i][v][g]=0;
                R[i1][v-1][g]=0;
                P[N_z[k].i][v-1][g]=1;
                R[N_z[k].i][v-1][g]=R[N_z[k].i][v][g-1];
                R[i1][v][g]=R[i1][v-1][g-1];

                v=v-1;

                log << "Perestanovka: \n";
                log << "s:" << s << " g:" << g << std::endl;
                log << "P: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        log << P[i][j][g] << " ";
                    }
                    log << std::endl;
                }
                log << "R: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        log << R[i][j][g] << " ";
                    }
                    log << std::endl;
                }
                step=8;
                break;
            }
            case 8:
            {
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        P_temp[i][j]=P[i][j][g];
                        R_temp[i][j]=R[i][j][g];
                    }
                }
                //вычисление матрицы t_jq_0l(s+g)
                start_time(P_temp, R_temp, n_z, n_p_z);
                //вычисление матрицы t_ji_nl(s+g)
                start_part(P_temp, R_temp, n_z, n_p_z);
                //вычисление критерия f3(s+g)
                f3[g]=objective(P_temp, R_temp, n_z, n_p_z);
                log << std::endl;
                log << "Kriterii (f3): " << f3[g] << std::endl;
                step=9;
                break;
            }
            case 9:
            {
                log << std::endl;
                log << "Effektivnoe zna4enie kriteriya (f3_eff): " << f3_eff << std::endl;
                log << "Novoe zna4enie kriteriya (f3): " << f3[g] << std::endl;

                if (f3_eff > f3[g])
                {
                    log << "Kriterii ulu4shilsya! \n";
                    //текущее решение более эффективное
                    for (int i=1; i<=n_z; i++)
                    {
                        for (int j=1; j<=n_p_z; j++)
                        {
                            P_eff[i][j]=P[i][j][0]=P[i][j][g];
                            R_eff[i][j]=R[i][j][0]=R[i][j][g];
                        }
                    }
                    f3_eff=f3[0]=f3[g];
                    //фиксируется локально-эффективное решение

                    log << "--------Lokalno-effektivnoe reshenie: s=" << s;
                    ofstream P_file(route_P);
                    log << " g=" << g << std::endl << "P: " << std::endl;
                    for (int i=1; i<=n_z; i++)
                    {
                        for (int j=1; j<=n_p_z; j++)
                        {
                            log << P_eff[i][j] << " ";
                            P_file << P_eff[i][j] << "\t";
                        }
                        log << std::endl;
                        P_file << "\n";
                    }
                    log << std::endl;
                    ofstream R_file (route_R);
                    log << "R: " << std::endl;
                    for (int i=1; i<=n_z; i++)
                    {
                        for (int j=1; j<=n_p_z; j++)
                        {
                            log << R_eff[i][j] << " ";
                            R_file << R_eff[i][j] << "\t";
                        }
                        log << std::endl;
                        R_file << "\n";
                    }
                    log << std::endl;
                    s++;
                    g=1;
                    if (v>1)
                    {
                        step=7;
                    }
                    else
                    {
                        step=13;
                    }
                }
                else
                {
                    step=10;
                }
                break;
            }
            case 10:
            {
                //проверка метрики
                w[g]=0;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1;j<=n_p_z2; j++)
                    {
                        w[g]+=abs(P[i][j][g]-P[i][j][0]);
                    }
                }
                w[g]/=2;
                log << std::endl << "w=" << w[g] << std::endl;
                if (w[g]<=w_max)
                {
                    //дальнейший поиск более эффективного решения
                    step=11;
                }
                else
                {
                    //получено локально-эффективное решение
                    step=12;
                }

                break;
            }
            case 11:
            {
                g++;
                if (v>1)
                {
                    step=7;
                }
                else
                {
                    step=13;
                }
                break;
            }
            case 12:
            {
                //фиксируется локально-эффективное решение

                log << "--------Lokalno-effektivnoe reshenie: s=" << s;
                ofstream P_file(route_P);
                log << " g=" << g << std::endl << "P: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        log << P_eff[i][j] << " ";
                        P_file << P_eff[i][j] << "\t";
                    }
                    log << std::endl;
                    P_file << "\n";
                }
                log << std::endl;
                ofstream R_file (route_R);
                log << "R: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        log << R_eff[i][j] << " ";
                        R_file << R_eff[i][j] << "\t";
                    }
                    log << std::endl;
                    R_file << "\n";
                }
                log << std::endl;
                step=13;
                break;
            }
            case 13:
            {
                mi--;
                h++;
                n_p_z2++;
                if (mi>0)
                {
                    s++;
                    v_max++;
                    step=2;
                }
                else
                {
                    step=14;
                }
                break;
            }
            case 14:
            {
                k++;
                if (k<=k_z)
                {
                    step=1;
                }
                else
                {
                    step=15;
                }
                break;
            }
            case 15:
            {
                stop=1;
                break;
            }
        }

    }

    if (k_z == 1)
    {
        //фиксируется локально-эффективное решение
        log << "--------Lokalno-effektivnoe reshenie: s=" << s;
        ofstream P_file(route_P);
        log << " g=" << g << std::endl << "P: " << std::endl;
        for (int i=1; i<=n_z; i++)
        {
            for (int j=1; j<=n_p_z; j++)
            {
                log << P_eff[i][j] << " ";
                P_file << P_eff[i][j] << "\t";
            }
            log << std::endl;
            P_file << "\n";
        }
        log << std::endl;
        ofstream R_file (route_R);
        log << "R: " << std::endl;
        for (int i=1; i<=n_z; i++)
        {
            for (int j=1; j<=n_p_z; j++)
            {
                log << R_eff[i][j] << " ";
                R_file << R_eff[i][j] << "\t";
            }
            log << std::endl;
            R_file << "\n";
        }
        log << std::endl;
    }

    log.close();
    iteration++;
    //Очистка памяти
  /*  for (int i=1; i<=k_z;i++)
    {
        delete (N_z[i].Ai);
        N_z[i].Ai=NULL;
    }*/
    delete [] N_z;
    delete[] P;
    delete[] R;
    delete[] P_eff;
    delete[] R_eff;
    delete[] R_temp;
    delete[] P_temp;
    delete[] f3;
    P=R=NULL;
    P_eff=R_eff=P_temp=R_temp=NULL;
    f3=NULL;
    //return a.exec();
    return 0;
}
