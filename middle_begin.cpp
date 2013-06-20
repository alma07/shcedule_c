#include <QCoreApplication>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string.h>


//ИНДЕКСАЦИЯ ВСЕХ МАССИВОВ С 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/*константы*/
const int II=4;          //количество типов
const int L=2;          //количество приборов
const int QQ=10;         //максимальное количество требований в партии
const char* route_til= "С:\\input\\time\\til.txt";
const char* route_tikl="С:\\input\\time\\tikl.txt";
const char* route_tiil="С:\\input\\time\\tiil.txt";
const char* route_tjq0l="С:\\input\\time\\tjq0l.txt";
const char* route_tjinl="С:\\input\\time\\tjinl.txt";
const char* route_middle="С:\\input\\middle.txt";
const char* route_P="С:\\input\\P.txt";
const char* route_R="С:\\input\\R.txt";

/*матрицы времён*/
float t_i_l[II+1][L+1],          //обработки требования i-го типа на l-ом приборе
      t_ik_l[II+1][II+1][L+1],    //переналадки l-го приб. с обр-ки i-го типа на k-ый
      t_ii_l[II+1][L+1];         //первонач. наладки l-го приб. на обраб i-го типа


//-----------------------------Функции---------------------------
void start_time (int **P, int **R, int n_z, int n_p_z);
float objective (int **P, int **R, int n_z, int n_p_z);
void start_part (int **P, int **R, int n_z, int n_p_z);

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
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
    const int razmer=10;
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


    //----------------Тестовые входные данные----------------------
    //файл должен быть записан правильно!!!!
    k_z=0;
    std::ifstream in (route_middle);
    char *str = new char [1024];
    while (!in.eof())
    {
        in.getline(str,1024,'\n');
        k_z++;
    }
    in.close();
    in.open (route_middle);
    I_z=new int[k_z+1];
    N_z=new Group[k_z+1];
    k=0;
    while (k<k_z)
    {
        k++;
        in >> I_z[k];
        N_z[k].i=I_z[k];
        std::cout << "i:" << I_z[k];
        in >> N_z[k].mi;
        std::cout << " m:" << N_z[k].mi;
        N_z[k].Ai=new int[N_z[k].mi];
        for (int i=1; i<=N_z[k].mi; i++)
        {
            in >> N_z[k].Ai[i];
        }
        std::cout << " A:";
        for (int i=1; i<=N_z[k].mi;i++)
        {
            std::cout << " " << N_z[k].Ai[i];
        }
        std::cout << std::endl;
    }
    in.close();
    std::cout << "Number of types in group: " << k << std::endl;

    //---------------Входные данные (матрицы времён)---------------
    std::ifstream til (route_til);
    for (int i=1; i<=II; i++)
    {
        for (int l=1; l<=L; l++)
        {
            til >> t_i_l[i][l];
        }
    }
    til.close();

    char c;
    std::ifstream tikl (route_tikl);
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

    std::ifstream tiil (route_tiil);
    for (int i=1; i<=II; i++)
    {
        for (int l=1; l<=L; l++)
        {
            tiil >> t_ii_l[i][l];
        }
    }
    tiil.close();


    //------------Начальная инициализация параметров-----------------
    g_max=20;
    k=1;
    v_max=0;
    n_z=k_z;

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

    std::cout << std::endl << "P: " << std::endl;
    for (int i=1; i<=n_z; i++)
    {
        for (int j=1; j<=n_p_z; j++)
        {
            std::cout << P_eff[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << "R: " << std::endl;
    for (int i=1; i<=n_z; i++)
    {
        for (int j=1; j<=n_p_z; j++)
        {
            std::cout << R_eff[i][j] << " ";
        }
        std::cout << std::endl;
    }

    k=2;
    while (stop!=1)
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
                std::cout << std::endl;
                std::cout << "P:" << std::endl;
                for (int i=1;i<=n_z;i++)
                {
                    for (int j=1; j<=n_p_z;j++)
                    {
                        std::cout << P[i][j][g] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
                std::cout << "R:" << std::endl;
                for (int i=1;i<=n_z;i++)
                {
                    for (int j=1; j<=n_p_z;j++)
                    {
                        std::cout << R[i][j][g]<< " ";
                    }
                    std::cout << std::endl;
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
                std::cout << std::endl;
                std::cout << "function: " << f3[g] << std::endl;

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

                std::cout << "g: " << g << std::endl;
                std::cout << "P: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        std::cout << P[i][j][g] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "R: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        std::cout << R[i][j][g] << " ";
                    }
                    std::cout << std::endl;
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
                std::cout << std::endl;
                std::cout << "function: " << f3[g] << std::endl;
                step=9;
                break;
            }
            case 9:
            {
                std::cout << std::endl;
                std::cout << "f3_eff: " << f3_eff << std::endl;
                std::cout << "f3[g]: " << f3[g] << std::endl;

                if (f3_eff > f3[g])
                {
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
                std::cout << std::endl << "w=" << w[g] << std::endl;
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

                std::cout << "--------Lokalno-effektivnoe reshenie: s=" << s;
                std::ofstream P_file(route_P);
                std::cout << " g=" << g << std::endl << "P: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        std::cout << P_eff[i][j] << " ";
                        P_file << P_eff[i][j] << "\t";
                    }
                    std::cout << std::endl;
                    P_file << "\n";
                }
                std::cout << std::endl;
                std::ofstream R_file (route_R);
                std::cout << "R: " << std::endl;
                for (int i=1; i<=n_z; i++)
                {
                    for (int j=1; j<=n_p_z; j++)
                    {
                        std::cout << R_eff[i][j] << " ";
                        R_file << R_eff[i][j] << "\t";
                    }
                    std::cout << std::endl;
                    R_file << "\n";
                }
                std::cout << std::endl;
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

    //Очистка памяти
    for (int i=1; i<=k_z;i++)
    {
        delete (N_z[i].Ai);
        N_z[i].Ai=NULL;
    }
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
    return a.exec();
}

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
                        t_jq_0l[j][1][l]=std::max(sum1,sum2);
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
                         t_jq_0l[j][q][l]=std::max(sum1,sum2);
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
            for (q=1; q<=Q; q++)
            {
                if (t_jq_0l[j][q][l]>0)
                {
                    tjq0l << t_jq_0l[j][q][l] << "\t";
                    std::cout << t_jq_0l[j][q][l] << " ";
                }
                else
                {
                    tjq0l << "0\t";
                    std::cout << "0 ";
                }
            }
            tjq0l << "\n";
            std::cout << std::endl;
        }
        tjq0l << ";\n";
        std::cout << ";" << std::endl;
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
            for (q=1; q<=Q; q++)
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


