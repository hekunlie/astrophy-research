#include<FQlib.h>

void task_alloc_d(const int total_task_num, const int division_num, const int my_part_id, int &my_st_id, int &my_ed_id, int *task_count, int *entry_for_gather)
{
    int i,j,m,n;
    m = total_task_num/division_num;
    n = total_task_num%division_num;

    for(i=0;i<division_num;i++)
    {
        task_count[i] = m;
        if(i<n){task_count[i] +=1;}
    }
    m=0;
    n=0;
    for(i=0;i<division_num;i++)
    {   
        m=0;
        for(j=0;j<i;j++)
        {
            m+=task_count[j];
        }
        entry_for_gather[i]=m;
    }
    n = m+task_count[my_part_id];
    my_st_id = entry_for_gather[my_part_id];
    my_ed_id = my_st_id + task_count[my_part_id];
}

int main(int argc, char **argv)
{
    int i,j,k;
    int numprocs, task_num;
    int my_st, my_ed;
    int *task_count, *entry_for_gather;

    numprocs = atoi(argv[1]);
    task_num = atoi(argv[2]);

    task_count = new int[numprocs];
    entry_for_gather = new int[numprocs];


    for(i=0;i<numprocs;i++)
    {   
        task_alloc(task_num, numprocs, i, my_st, my_ed, task_count,entry_for_gather);
        std::cout<<i<<" "<<my_st<<" "<<my_ed<<std::endl;
        if(i==0){show_arr(task_count,1,numprocs); show_arr(entry_for_gather,1,numprocs);}
        std::cout<<std::endl;
    }

    return 0;
}