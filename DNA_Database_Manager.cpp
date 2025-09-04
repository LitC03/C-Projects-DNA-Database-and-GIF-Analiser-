/*
Assignment 2 -  DNA Editing
By Lito Paraskevi Chatzidavari
CID: 01711772
26-03-2021

This program takes .fna files and stores them into a database. The user can then use the files saved to find
other sequences either by input or other FASTA files. Modifications can also be made and are inmediately done to the database.
However, if a sequence were to be deleted completely, no further modifications will be available.
The program can also analyse all the database to search for a specific sequence.
*/

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include "tictoc.h"

using namespace std;

//-----------------------------------------------------------------------------

struct string_for_class
{
  char* p;
  int size;
};

//-----------------------------------------------------------------------------

class Nucleotide //Nucleotide is equivalent to a "node" of the Index linked list
{
  private:
    char n;
    Nucleotide* n_next;
  public:
    Nucleotide();
    Nucleotide(char n_in, Nucleotide* n_next_in);
    char getN();
    void setN(char n_in);
    Nucleotide* getNextN();
    void setNextN(Nucleotide* n_next_in);

};
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

class Index //Index is equivalent to a "node" of the DNADatabase linked list
{
  private:
    int x;
    Index* x_next;
    Nucleotide* n_head;
    Nucleotide* n_last; //pointer to the last nucleotide of the sequence
    string filename; //file from which information was extracted
    int nucl_num; //total number of nucleotides in Index
  public:
    Index();
    Index(int x_in, Index* x_next_in, string filename_in);
    int getX();
    string get_filename();
    Index* getNext();
    int get_nucl_num();
    Nucleotide* get_n_head();
    void setNext(Index* x_next_in);
    void setX(int x_in);
    void setname(string filename_in);
    void addhere(Nucleotide* here, char c, int position, int length);
    void delete_here(Nucleotide* previous, int position, int length);
    void popN();
    friend class DNADatabase;
    ~Index();
};
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

class DNADatabase
{
  private:
    Index* p_head;
    int size;
    Index* p_last;//pointer to the last Index of the database
  public:
    DNADatabase();
    int getsize();
    void setsize(int size_in);
    Index* get_p_head();
    void push(string filename_in);
    void push(Index* p_new);
    void addbegin(string filename_in);
    void find(string_for_class sequence, string file_to_search);
    void pop_begin();
    void print();
    int find_from_file(string file_to_open, string file_to_search, const bool& analysis, int matches);
    void add_NUCL(string_for_class sequence, int position, string file_to_search, bool& changes_done);
    void add_NUCL_from_file(string file_to_add, int position, string filename, bool& changes_done);
    void deleteN(int position, int length, string file_to_search, bool& changes_done);
    void change(int position, int length, string_for_class sequence, string file_to_search, bool& changes_done);
    void change_from_file(int position, int length, string filename, string file_to_search, bool& changes_done);
    void save(string file_to_save, string filename);
    void analysis(string file_find, const vector<string>& filename_list);
    ~DNADatabase();
};
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void option1(vector<string>& filename_list, DNADatabase& dna_db);
string get_filename(int position, string& reverse);
bool check_file(string& filename);
void option2(vector<string>& filename_list, DNADatabase& dna_db);
void option3(vector<string>& filename_list, DNADatabase& dna_db);
bool check_string_is_num(string& input_str);
void submenu_opt2(const string& filename, DNADatabase& dna_db);
void option2_1(const string& filename, DNADatabase& dna_db);
void option2_2(const string& filename, DNADatabase& dna_db);
void option2_3(const string& filename, DNADatabase& dna_db, bool& changes_done);
void option2_4(const string& filename, DNADatabase& dna_db, bool& changes_done);
void option2_5(const string& filename, DNADatabase& dna_db, bool& changes_done);
void option2_6(const string& filename, DNADatabase& dna_db, bool& changes_done);
void option2_7(const string& filename, DNADatabase& dna_db, bool& changes_done);
void option2_8(DNADatabase& dna_db, string filename, const bool& changes_done);
string_for_class string_to_chars(string str);
void print_find(Index* p_seek, Nucleotide* n_compare, int i, int k, int nucl_read, Nucleotide* n_seek, const bool& analysis);
void print_before(Nucleotide* n_seek, int position);
void print_after(Nucleotide* n_seek, int position, int nucl_number);
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

Nucleotide::Nucleotide()
{
  n = '\0';
  n_next = nullptr;
  cout << "An empty nucleotide has been added" << endl;
}

Nucleotide::Nucleotide(char n_in, Nucleotide* n_next_in)
{
  n = n_in;
  n_next = n_next_in;
}

char Nucleotide::getN()
{
  return n;
}

void Nucleotide::setN(char n_in)
{
  n = n_in;
}


Nucleotide* Nucleotide::getNextN()
{
  return n_next;
}


void Nucleotide::setNextN(Nucleotide* n_next_in)
{
  n_next = n_next_in;
}

//-----------------------------------------------------------------------------
Index::Index()
{
  x = 0;
  x_next = nullptr;
  n_head = nullptr;
  n_last = nullptr;
  filename = "";
  nucl_num = 0;
}

Index::Index(int x_in, Index* x_next_in, string filename_in)
{
  fstream file;
  char c;
  Nucleotide *n_new;

  nucl_num = 0;
  x = x_in;
  x_next = x_next_in;
  filename = filename_in;
  n_head = nullptr;

  //File is opened and every character after the first line is put into the Index we are working on.
  file.open(filename, fstream::out | fstream::in);

  while(c!='\n')
  {
    c = file.get();
  }

  while(c!=EOF)
  {
    c = file.get();
    if ((c!='\n')&&(c!='\0')&&(c!=EOF))
    {
      //the number of nucleotides in the sequence increases with every new nucleotide
      nucl_num++;
      n_new = new Nucleotide(c, nullptr);

      //every new nucleotide is put at the last position of the Index
      if (n_head == nullptr)
        n_head = n_last = n_new;
      else
      {
        n_last->setNextN(n_new);
        n_last=n_last->getNextN();
      }
    }
  }
}

void Index::addhere(Nucleotide* here, char c, int position, int length)
{
  //the pointer "here" is found 2 positions before where the new nucleotide will be added
  Nucleotide *n_new , *next_n, *n_seek;

  if (here==nullptr)
  {
    // add nucleotide at the last position

    n_new = new Nucleotide(c, nullptr);
    n_last->setNextN(n_new);
    n_last=n_last->getNextN();

  }
  //there is a special case when the user wants to add nucleotides at the beggining of the sequence
  else if(position==0)
  {
    if (length==0)
    {
      //If we are adding only one nucleotide (or the first of the nucleotides we want to add), we can treat this situation as a queue
      n_head = new Nucleotide(c, n_head);
    }
    else
    {
      //After a new nucleotide has been added, the whole sequence has been shifted. This means that "here" is no longer of use.
      // We need a new pointer, which is found using the number of nucleotides that have already been added
      n_seek = n_head;
      for(int i=0; i< length-1; i++)
         n_seek = n_seek->getNextN();

      next_n = n_seek->getNextN();
      n_new = new Nucleotide(c, next_n);
      n_seek->setNextN(n_new);
    }
  }
  else
  {
    //add nucleotide somewhere in the middle of the index
    next_n = here->getNextN();

    n_new = new Nucleotide(c, next_n);
    here->setNextN(n_new);
  }
  nucl_num++;
}

void Index::delete_here(Nucleotide* previous, int position, int length)
{
  //the pointer "previous" is found in the nuleotide 2 positions before the nucleotide that needs to be deleted
  //the pointer "here" points to the nucleotide that will be deleted
  //the pointer "after_del" is the one found inside the nucleotide that will be deleted
  Nucleotide* here = previous->getNextN();
  Nucleotide* after_del;
  Nucleotide* n_seek;

  if (nucl_num==0)
    cout << "Index is empty." << endl;
  else if (here==nullptr)
    cout << "this is not possible" << endl;
  else
  {
    //Again there is a special case if we want to delete the first nuleotide of the sequence
    if (position==0)
    {
      //If we are deleting only one nucleotide (or the first of the nucleotides we want to delete), we just need to change n_head
      if (length==0)
        n_head=here;
      else
      {
        //If have arleady deleted 1 nucleotide, n_head needs to point to the second nucleotide of the unmodified sequence
        n_head = n_head->getNextN();
      }
    }
    else
    {
      after_del = here->getNextN();

      //if the nucleotide being deleted is the last one, n_last has to be changed
      if (after_del==nullptr)
      {
        n_last = previous;
        previous->setNextN(nullptr);
      }
      else
      {
        previous->setNextN(after_del);
        here->setNextN(nullptr);
      }
    }
    //the number of nucleotides in the sequence decreases
    nucl_num--;
  }
}

Index* Index::getNext()
{
  return x_next;
}

int Index::getX()
{
  return x;
}

string Index::get_filename()
{
  return filename;
}

int Index::get_nucl_num()
{
  return nucl_num;
}

Nucleotide* Index::get_n_head()
{
  return n_head;
}

void Index::setNext(Index* x_next_in)
{
    x_next = x_next_in;
}

void Index::setX(int x_in)
{
  x = x_in;
}

void Index::setname(string filename_in)
{
  filename = filename_in;
}

void Index::popN()
{
  Nucleotide* n_temp;
  n_temp = n_head->getNextN();
  delete n_head;
  n_head = n_temp;
  nucl_num--;
}

Index::~Index()
{
  while (n_head != nullptr)
  {
    popN();
  }
}

//-----------------------------------------------------------------------------

DNADatabase::DNADatabase()
{
  p_head = nullptr;
  size = 0;
  p_last = nullptr;
}

int DNADatabase::getsize()
{
  return size;
}

Index* DNADatabase::get_p_head()
{
  return p_head;
}

void DNADatabase::setsize(int size_in)
{
  size = size_in;
}

void DNADatabase::addbegin(string filename_in)
{
  Index* p_seek;
  //new index is added at the beggining of the database
  p_head = new Index(1, p_head, filename_in);
  size++;
  p_seek = p_head;
  //enumeration is shifted
  for (int i=1; i<size+1; i++)
  {
    p_seek->setX(i);
    p_seek=p_seek->getNext();
  }
}

void DNADatabase::push(string filename_in)
{
  size++;
  Index* p_new;

  p_new = new Index(size, nullptr, filename_in);

  //if database is empty, p_head needs to be changed when a sequence is added
  if (p_head == nullptr)
    p_head = p_last = p_new;
  else
  {
    //if not, an index is added at the end of the database
    p_last->setNext(p_new);
    p_last=p_last->getNext();
  }
}

void DNADatabase::push(Index* p_new)
{
  size++;

  if (p_head == nullptr)
    p_head = p_last = p_new;
  else
  {
    p_last->setNext(p_new);
    p_last=p_last->getNext();
  }
}


void DNADatabase::pop_begin()
{
  Index* p_temp;

  p_temp = p_head->getNext();
  delete p_head;
  p_head = p_temp;

  //for the sake of completeness, this should uncommented if it were to be used somewhere in the code (and not just the destructor)
  // for (int i=1; i<size; i++)
  // {
  //   p_seek->setX(i);
  //   p_seek=p_seek->getNext();
  // }
  size--;
}

void DNADatabase::print()
{
  Index* p_seek = p_head;
  Nucleotide* n_seek;
  int x_out;
  int nucl_out;
  char nucl;
  string file;

  cout << "\nThe elements in the list are:" << endl;
  for (int i=0; i<size; i++)
  {
    //get all information of an index
    x_out = p_seek->getX();
    file = p_seek->get_filename();
    nucl_out = p_seek->get_nucl_num();
    n_seek = p_seek->get_n_head();

    cout << "\n" << x_out << ". " << file << endl;

    cout << "It contains " << nucl_out << " bases, which are:" << endl;

    for (int j=0; j<nucl_out; j++)
    {
      nucl = n_seek->getN();
      cout << nucl;
      n_seek = n_seek->getNextN();
    }

    cout << endl;

    p_seek = p_seek->getNext();
  }
}

void DNADatabase::find(string_for_class sequence, string file_to_search)
{
  Index* p_seek;
  int nucl_read;
  int matches=0;
  bool comparison;
  Nucleotide *n_seek, *n_compare;

  for(p_seek = p_head; file_to_search!=p_seek->get_filename(); p_seek = p_seek->getNext());

  n_seek = p_seek->get_n_head();
  nucl_read = p_seek->get_nucl_num();

  //if only one character needs to be found, every time the character is found is considered a match
  if (sequence.size==1)
  {
    for(int i=0; i<nucl_read; i++)
    {
      if (n_seek->getN()==sequence.p[0])
      {
        matches++;
        cout << "\nMatch #" << matches-1 << endl;
        print_find(p_seek, n_seek, i, 0, nucl_read , n_seek, false);
      }
      n_seek=n_seek->getNextN();
    }
  }
  else
  {
    for(int i=0; i<nucl_read-1; i++)
    {
      //If the nucleotide being read equals the first element of the sequence, we check if the next nucleotides match as well
      if (n_seek->getN()==sequence.p[0])
        {
          comparison = true;
          n_compare = n_seek->getNextN();
          for(int k=1; ((k<sequence.size)&&(comparison)); k++)
          {
            if(n_compare->getN()!=sequence.p[k])
            {
              comparison = false;
            }
            //if all nucleotides match, the result is printed
            else if ((k==sequence.size-1)&&(n_compare->getN()==sequence.p[sequence.size-1])&&(comparison))
              {
                matches++;
                cout << "\nMatch #" << matches-1 << endl;
                print_find(p_seek, n_compare, i, k, nucl_read , n_seek, false);
              }
            if (n_compare->getNextN()!=nullptr)
              n_compare = n_compare->getNextN();
            else comparison = false;
          }
        }
      n_seek = n_seek->getNextN();
    }
  }
}

int DNADatabase::find_from_file(string file_to_open, string file_to_search, const bool& analysis, int matches_an)
{
  char c, from_file, searching;
  Index* file_p_head;
  Index* p_seek;
  bool comparison;
  int nucl_read, nucl_opened_file;
  int matches=0;
  Nucleotide *n_opened_file, *n_seek, *n_compare;


  for(p_seek = p_head; file_to_search!=p_seek->get_filename(); p_seek = p_seek->getNext());

  n_seek = p_seek->get_n_head();
  nucl_read = p_seek->get_nucl_num();

  //a new Index object is made from the file
  file_p_head = new Index(1, nullptr, file_to_open);

  nucl_opened_file = file_p_head->get_nucl_num();


  if(nucl_opened_file!=0)
  {
    //same algorithm as before is followed, this time, nucleotides are being read from an index object
    if (nucl_opened_file==1)
    {
      for(int i=0; i<nucl_read; i++)
      {
        n_opened_file = file_p_head->get_n_head();
        if (n_seek->getN()==n_opened_file->getN())
        {
          matches++;
          cout << "\nMatch #";
          if (analysis)
          {
            cout << matches_an << endl;
            matches_an++;
          }
          else
            cout << matches-1 << endl;

          print_find(p_seek, n_seek, i, 0, nucl_read , n_seek, analysis);
        }
        n_seek=n_seek->getNextN();
      }
    }
    else
    {
      for(int i=0; i<nucl_read-1; i++)
      {
        n_opened_file = file_p_head->get_n_head();
        //If the nucleotide being read equals the first element of the sequence, we check if the next nucleotides match as well
        if (n_seek->getN()==n_opened_file->getN())
          {
            comparison = true;
            n_compare = n_seek->getNextN();

            for(int k=1; ((k<nucl_opened_file)&&(comparison)); k++)
            {
              n_opened_file = n_opened_file->getNextN();
              if(n_compare->getN()!=n_opened_file->getN())
              {
                comparison = false;
              }
              //if a match is found the results are printed
              else if ((k==nucl_opened_file-1)&&(n_compare->getN()==n_opened_file->getN())&&(comparison))
                {
                  matches++;
                  cout << "\nMatch #";
                  if (analysis)
                  {
                    cout << matches_an << endl;
                    matches_an++;
                  }
                  else
                    cout << matches-1 << endl;
                  print_find(p_seek, n_compare, i, k, nucl_read , n_seek, analysis);
                }
              if (n_compare->getNextN()!=nullptr)
                n_compare = n_compare->getNextN();
              else comparison = false;
            }
          }
        n_seek = n_seek->getNextN();
      }
    }
  }
  return matches_an;
}

void DNADatabase::add_NUCL(string_for_class sequence, int position, string file_to_search, bool& changes_done)
{
  Index *p_seek;
  int nucl_number;
  fstream file;

  //find pointer to corresponding Index
  for(p_seek = p_head; file_to_search!=p_seek->get_filename(); p_seek = p_seek->getNext());

  nucl_number = p_seek->get_nucl_num();

  if(nucl_number==0)
  {
    cout << "\nIt appears that your sequence is empty. Please, try with a new one." << endl;
  }
  else if (position>nucl_number)
  {
    cout << "\nThe position you have selected is unavailable. Please try again." << endl;
  }
  else
  {
    //a file with the sequence that will be added is made and it is passed to the function for adding nucleotides from a file
    file.open("DNA.txt",fstream::out);
    if (file.is_open())
    {
      file.put('>');
      file.put('\n');
      for (int i=0; i<sequence.size; i++)
        file.put(sequence.p[i]);
    }
    else cout << "There was an error" << endl;

    file.close();

    add_NUCL_from_file("DNA.txt", position, file_to_search, changes_done);

    //the file is removed from the folder
    remove("DNA.txt");
  }
}

void DNADatabase::add_NUCL_from_file(string file_to_add, int position, string filename, bool& changes_done)
{
  Index *new_seq_head, *p_seek;
  Nucleotide* n_seek, *n_adding , *n_help, *n_compare;
  int nucl_number;

  //find pointer to corresponding Index
  for(p_seek = p_head; filename!=p_seek->get_filename(); p_seek = p_seek->getNext());

  nucl_number = p_seek->get_nucl_num();

  if(nucl_number==0)
  {
    cout << "\nIt appears that your sequence is empty. Please, try with a new one." << endl;
  }
  else if (position>nucl_number)
  {
    cout << "\nThe position you have selected is unavailable. Please try again." << endl;
  }
  else
  {
    //new object index is made from file
    new_seq_head = new Index(1, nullptr,file_to_add);

    n_seek = p_seek->get_n_head();
    n_adding = new_seq_head->get_n_head();

    //get pointer to where nucleotides need to be added
    for(int i=0; i<position-1; i++)
    {
      n_seek = n_seek->getNextN();
    }

    //add nucleotides as they are read from the new object made recently
    for (int i=0; i<new_seq_head->get_nucl_num(); i++)
    {
      p_seek->addhere(n_seek, n_adding->getN(), position, i);
      n_seek = n_seek->getNextN();
      n_adding = n_adding->getNextN();
    }

    //pointers are found to print results correctly
    n_help = p_seek->get_n_head();

    for(int i=0; i<position; i++)
      n_help = n_help->getNextN();

    n_compare = n_help;
    for(int i=0; i<new_seq_head->get_nucl_num()-1; i++)
      n_compare = n_compare->getNextN();

    //results are printed
    cout << endl;
    print_find(p_seek, n_compare ,position, new_seq_head->get_nucl_num()-1, p_seek->get_nucl_num(), n_help, false);
    changes_done = true;
  }
}

void DNADatabase::deleteN(int position, int length, string file_to_search, bool& changes_done)
{
  Index *p_seek, *not_saved_p_head;
  Nucleotide *n_seek, *n_help, *n_compare;
  int nucl_number;
  char nucl;

  //find pointer to corresponding Index
  for(p_seek = p_head; file_to_search!=p_seek->get_filename(); p_seek = p_seek->getNext());

  nucl_number = p_seek->get_nucl_num();

  if(nucl_number==0)
  {
    cout << "\nIt appears that your sequence is empty. Please, try with a new one." << endl;
  }
  else if(position>nucl_number)
  {
    cout << "\nThe position you have selected is unavailable. Please try again." << endl;
  }
  else if ((position+length>nucl_number)||((position==0)&&(length==nucl_number-1)))
  {
    cout << "\nYou are trying to delete unavailable nucleotides. Please try again." << endl;
  }
  else
  {
    //get pointers to print situation before deletion, and print it
    n_seek = p_seek->get_n_head();

    for(int i=0; i<position; i++)
    {
      n_seek = n_seek->getNextN();
    }
    n_compare = n_seek;
    for(int i=0; i<length-1; i++)
    {
      n_compare = n_compare->getNextN();
    }
    cout << "\nDNA sequence deletion information:" << endl;
    print_find(p_seek, n_compare, position, length-1, nucl_number, n_seek, false);

    //get pointer to where nucleotides need to be deleted
    n_seek = p_seek->get_n_head();

    for(int i=0; i<position-1; i++)
    {
      n_seek = n_seek->getNextN();
    }
    //delete nucleotides
    for (int i=0; i<length; i++)
    {
      p_seek->delete_here(n_seek, position, i);
    }

    cout << "\nThe DNA sequence has been deleted." << endl;

    cout << "\nDNA sequence deletion result:" << endl;
    cout << "base pair position:\t" << position << endl;

    //print 10 nucleotides before deleted sequence
    n_seek = p_seek->get_n_head();
    print_before(n_seek, position);

    for(int i=0; i<position; i++)
       n_seek = n_seek->getNextN();

    //print 10 nucleotides after deleted sequence
    print_after(n_seek, position+length, nucl_number);
    changes_done=true;
  }
}

void DNADatabase::change(int position, int length, string_for_class sequence, string file_to_search, bool& changes_done)
{
  Index *p_seek;
  int nucl_number;
  fstream file;

  //find pointer to corresponding Index
  for(p_seek = p_head; file_to_search!=p_seek->get_filename(); p_seek = p_seek->getNext());

  nucl_number = p_seek->get_nucl_num();

  if(nucl_number==0)
  {
    cout << "\nIt appears that your sequence is empty. Please, try with a new one." << endl;
  }
  else if (position>nucl_number)
  {
    cout << "\nThe position you have selected is unavailable. Please try again." << endl;
  }
  else if (position+length>nucl_number)
  {
    cout << "\nYou are trying to delete unavailable nucleotides. Please try again." << endl;
  }
  else
  {
    //new file is created with the sequence entered by the user
    file.open("DNA.txt",fstream::out);
    if (file.is_open())
    {
      file.put('>');
      file.put('\n');
      for (int i=0; i<sequence.size; i++)
      {
        file.put(sequence.p[i]);
      }
    }
    else cout << "There was an error" << endl;

    file.close();

    //changes are made with "change from file" function
    change_from_file(position, length, "DNA.txt", file_to_search, changes_done);

    remove("DNA.txt");
  }
}

void DNADatabase::change_from_file(int position, int length, string filename, string file_to_search, bool& changes_done)
{
  Index *p_seek, *new_seq_head;
  Nucleotide *n_seek, *n_help, *n_compare, *n_adding;
  int nucl_number;
  fstream file;

  //find pointer to corresponding Index
  for(p_seek = p_head; file_to_search!=p_seek->get_filename(); p_seek = p_seek->getNext());

  nucl_number = p_seek->get_nucl_num();

  if(nucl_number==0)
  {
    cout << "\nIt appears that your sequence is empty. Please, try with a new one." << endl;
  }
  else if (position>nucl_number)
  {
    cout << "\nThe position you have selected is unavailable. Please try again." << endl;
  }
  else if (position+length>nucl_number)
  {
    cout << "\nYou are trying to delete unavailable nucleotides. Please try again." << endl;

  }
  else
  {
    //new Index object is created from file
    new_seq_head = new Index(1, nullptr,filename);

    n_seek = p_seek->get_n_head();
    n_adding = new_seq_head->get_n_head();

    //pointers are set to delete the replaced nucleotides
    for(int i=0; i<position-1; i++)
      n_seek = n_seek->getNextN();

    n_compare = n_seek;

    //nucleotides are deleted
    for (int i=0; i<length; i++)
    {
      p_seek->delete_here(n_seek, position, i);
    }

    //new nucleotides are added
    for (int i=0; i<new_seq_head->get_nucl_num(); i++)
    {
      p_seek->addhere(n_compare, n_adding->getN(), position, i);
      n_compare = n_compare->getNextN();
      n_adding = n_adding->getNextN();
    }

    if(position==0)
    {
      //get pointers to print results
      n_seek=p_seek->get_n_head();
      n_help = p_seek->get_n_head();

      for(int i=0; i<position; i++)
        n_help = n_help->getNextN();

      n_compare = n_help;
      for(int i=0; i<new_seq_head->get_nucl_num()-1; i++)
        n_compare = n_compare->getNextN();
    }
    else
      n_seek = n_seek->getNextN();

    //results are printed
    cout << endl;
    print_find(p_seek, n_compare, position, new_seq_head->get_nucl_num()-1, p_seek->get_nucl_num(), n_seek, false);
    changes_done=true;
  }
}

void DNADatabase::save(string file_to_save, string filename)
{
  fstream file;
  Nucleotide* n_seek;
  int count=0;
  int nucl_number;
  Index* p_seek;

  //find pointer to corresponding Index
  for(p_seek = p_head; file_to_save!=p_seek->get_filename(); p_seek = p_seek->getNext());

  n_seek = p_seek->get_n_head();
  nucl_number = p_seek->get_nucl_num();

  //new file is made with the characters contained in the edited Index
  file.open(filename,fstream::out);
  if (file.is_open())
  {
    file.put('>');
    file.put('\n');
    for (int i=0; i<nucl_number; i++)
    {
      file.put(n_seek->getN());
      n_seek = n_seek->getNextN();
      count++;
      // every 75 characters, a '\n' is added to make file more readable
      if (count==75)
      {
        file.put('\n');
        count=0;
      }
    }
    cout << "\nThe file has been added to your folder." << endl;
  }
  else cout << "\nThere was an error" << endl;

  file.close();
}

void DNADatabase::analysis(string file_find, const vector<string>& filename_list)
{
  int matches=0;

  //find_from_file is applied to every sequence of the database
  for (int i=0; i<size; i++)
  {
    matches = find_from_file(file_find, filename_list[i], true, matches);
  }
}

DNADatabase::~DNADatabase()
{
  while (p_head != nullptr)
  {
    pop_begin();
  }
}

//-----------------------------------------------------------------------------

int main()
{
  string input_str;
  vector<string> filename_list;
  DNADatabase dna_db;


  cout << "\nWelcome to the DNA Editing program" << endl;
  //display menu and ask for input
  while (input_str != "4")
  {
    cout << "\nSelect an option:" << endl;
    cout << "(1) Load DNA(s)." << endl;
    cout << "(2) Process a DNA." << endl;
    cout << "(3) Analyse the DNA database" << endl;
    cout << "(4) Quit." << endl;

    getline(cin,input_str);

    if (input_str == "1")
      option1(filename_list,dna_db);
    else if (input_str == "2")
      option2(filename_list, dna_db);
    else if (input_str == "3")
      option3(filename_list, dna_db);
    else if (input_str != "4")
      cout << "\n\"" << input_str << "\" is not a valid input. Please try again." << endl;
  }

  cout << "\nExiting program" << endl;

  return 0;
}

void option1(vector<string>& filename_list, DNADatabase& dna_db)
{
  string files_string;
  string again;
  int position, size;
  vector<string>::iterator it;
  string filename;
  bool valid_file;
  bool valid = true;
  TicToc time;

  while (valid)
  {
    position=-4;
    again = "1";

    cout << "\nLoad DNA strands" << endl;
    cout << "\nEnter the DNA file names:" << endl;
    cout << "For multiple files, separate them by a comma. Only .fna are recognised." << endl;
    cout << "Do not input anything and press enter if you wish to return to the menu." << endl;
    getline(cin,files_string);
    cout << endl;
    time.tic();
    size = files_string.length();

    if (files_string.empty())
    {
      valid = false;
      again = "2";
    }

    while (position < size)
    {
      //find position of ".fna" strings inside the user's input
      position = files_string.find(".fna",position+4);
      if (position != string::npos)
      {
        //if a match is found, get the filename
        filename = get_filename(position, files_string);
        if (!filename.empty())
        {
          filename += ".fna";
          // check if filename exists in folder and is in the FASTA format
          valid_file = check_file(filename);
          if (valid_file)
          {
            //find if file entered by user already exists in the program, if not, add it to database
            it = find (filename_list.begin(), filename_list.end(), filename);
            if (it != filename_list.end())

              cout << "\nThe file \"" << filename << "\" has already exists in the program. It won't be added again." << endl;
            else
            {
              cout << "Loading file #" << filename_list.size()+1 << " " << "\"" << filename << "\"..." << endl;
              dna_db.push(filename);
              filename_list.push_back(filename);
            }
          }
        }
      }
      else position = size;
    }

    time.toc();
    cout << endl;
    cout << time << endl;

    while (again != "2")
    {
      if (filename_list.empty())
        cout << "\nYou have not added any files since the start of the program." << endl;
      else
      {
        cout << "\nThe following files have been added since the start of the program:" << endl;
        for (int i=0; i<filename_list.size(); i++)
          cout << "#" << i+1 << " \"" << filename_list[i] << "\"" << endl;
      }

      cout << "\nWould you like to add more?" << endl;
      cout << "(1) Yes\n(2) No" << endl;
      getline(cin, again);
      if (again == "1")
        again = "2";
      else if (again == "2")
        valid = false;
      else if (again != "2")
        cout << "\n\"" << again << "\" is not a valid input. Please try again." << endl;
    }
  }
}

string get_filename(int position, string& files_string)
{
  //create empty string
  string filename("");
  string read;
  int spaces_num=0;
  string::iterator it_end = files_string.end();
  string::iterator it_check = files_string.begin()+position+4;

  //check if the last element of the string is the element after ".fna" or if ".fna" is followed by a coma or a space
  if ((it_check==it_end)||(files_string.substr(position+4, 1)==",")||(files_string.substr(position+4, 1)==" "))
  {
    //each character will be read from right to left until a comma is found or the string ends
    for(int i=position-1; (read!=",")&&(i>=0); i--)
    {
      read = files_string[i];
      //everything but commas are read
      if (read!=",")
      {
        //count spaces between characters
        if (read==" ")
          spaces_num++;
        else spaces_num = 0;

        //add characters to the string filename
        filename += read;
      }
    }
  }

  //delete preceding spaces
  if ((spaces_num)&&(!filename.empty()))
    filename.erase(filename.length()-spaces_num,spaces_num);

  //reverse filename since it has been read from right to left
  reverse(filename.begin(), filename.end());
  return filename;
}

bool check_file(string& filename)
{
  bool good;
  char c;
  fstream file;

  file.open(filename, fstream::in);

  //check if file exists in folder
  if (file.fail())
  {
    good = false;
    cout << "\nThere was an error opening the file \"" << filename << "\". It will not be used." << endl;
  }
  else
  {
    //check if file is in the FASTA format
    c = file.get();
    if (c == '>')
      good = true;
    else
    {
      good = false;
      cout << "\nThe file \"" << filename << "\" is not in the FASTA format. It will not be used." << endl;
    }
  }
  file.close();
  return good;
}

void option2(vector<string>& filename_list, DNADatabase& dna_db)
{
  bool show_menu2 = true;
  bool is_num;
  bool empty_list_message;
  string input_1st, file_to_process;
  string more;
  int file_index;

  cout << "\nProcess a DNA" << endl;

  while(show_menu2)
  {
    cout << "\nSelect a DNA to process:" << endl;
    if (filename_list.empty())
    {
      empty_list_message = true;
      while(empty_list_message)
      {
        cout << "\nThere are no available files at the moment." << endl;
        cout << "Do not input anything and press enter to return to the menu." << endl;
        getline(cin, input_1st);
        if (input_1st.empty())
        {
          empty_list_message = show_menu2 = false;
        }
        else
          cout << "\n\"" << input_1st << "\" is not a valid input. Please try again." << endl;
      }
    }
    else
    {
      for (int i=0; i<filename_list.size(); i++)
        cout  <<"(" << i+1 << ") " << filename_list[i] << endl;

      cout << "\nDo not input anything and press enter if you wish to return to the menu." << endl;
      getline(cin,input_1st);

      if (input_1st.empty())
        show_menu2 = false;
      else
      {
        is_num=check_string_is_num(input_1st);

        if (is_num)
        {
          // get number from users input
          file_index = stoi(input_1st);

          file_index--;

          if ((file_index < 0)||(file_index>=filename_list.size()))
            cout << "\n\"" << input_1st << "\" is out of range. Please try again" << endl;
          else
          {
            //if user has selected an available file, show submenu
            file_to_process = filename_list[file_index];
            cout << "\nYou have selected \"" << file_to_process << "\"" << endl;
            submenu_opt2(file_to_process, dna_db);

            more = "1";
            while (more != "2")
            {
              cout << "\nDo you want to process another DNA?" << endl;
              cout << "(1) Yes\n(2) No" << endl;
              getline(cin, more);
              if (more == "1")
                more = "2";
              else if (more == "2")
                show_menu2 = false;
              else if (more != "2")
                cout << "\n\"" << more << "\" is not a valid input. Please try again." << endl;
            }
          }
        }
        else
          cout << "\nRemember to enter numbers only." << endl;
      }
    }
  }
}

bool check_string_is_num(string& input_str)
{
  bool number = true;
  //check if there are spaces in the string
  int position = input_str.find(" ");

  //if there aren't any, check is all characters in it are digits or the number entered is a negative number
  if (position == string::npos)
    {
      if (all_of(input_str.begin(), input_str.end(), ::isdigit))
        number=true;
      else if ((all_of(input_str.begin()+1, input_str.end(), ::isdigit))&&(input_str.substr(0,1)=="-"))
        number=true;
      else
        number= false;
    }
  else
    number = false;

  if (input_str.empty())
    number = false;

  //if the number is too big for an int to contain, display error message
  if ((number)&&(input_str.length()>9))
  {
    cout << "\nThe number you have entered is too big in size. Please try again." << endl;
    number = false;
  }

  return number;
}

void submenu_opt2(const string& filename, DNADatabase& dna_db)
{
  string input_str = "1";
  bool changes_done=false;
  while (input_str!="9")
  {
    cout << "\nSelect from one of the following options" << endl;
    cout << "(1) Find DNA sequence by input" << endl;
    cout << "(2) Find DNA sequence by file" << endl;
    cout << "(3) Add DNA sequence by input" << endl;
    cout << "(4) Add DNA sequence by file" << endl;
    cout << "(5) Delete DNA sequence by input" << endl;
    cout << "(6) Replace DNA sequence by input" << endl;
    cout << "(7) Replace DNA sequence by file" << endl;
    cout << "(8) Save edited DNA sequence" << endl;
    cout << "(9) Exit submenu" << endl;
    getline(cin, input_str);

    if (input_str == "1")
      option2_1(filename, dna_db);
    else if (input_str == "2")
      option2_2(filename, dna_db);
    else if (input_str == "3")
      option2_3(filename, dna_db, changes_done);
    else if (input_str == "4")
      option2_4(filename, dna_db, changes_done);
    else if (input_str == "5")
      option2_5(filename, dna_db, changes_done);
    else if (input_str == "6")
      option2_6(filename, dna_db, changes_done);
    else if (input_str == "7")
      option2_7(filename, dna_db, changes_done);
    else if (input_str == "8")
      option2_8(dna_db, filename, changes_done);
    else if (input_str != "9")
        cout << "\n\"" << input_str << "\" is not a valid input. Please try again." << endl;
  }
}

void option2_1(const string& filename, DNADatabase& dna_db)
{
  string_for_class sequence;
  string user_input;
  TicToc time;

  cout << "\nFind DNA sequence by input" << endl;
  cout << "\nEnter the DNA sequence to search (eg, GTCACT):" << endl;
  getline(cin, user_input);

  time.tic();

  //get string to an array
  sequence = string_to_chars(user_input);

  dna_db.find(sequence, filename);
  time.toc();
  cout << endl;
  cout << time << endl;

  delete[] sequence.p;
}

string_for_class string_to_chars(string str)
{
  string_for_class a;
  int length = str.length();
  char* c = new char[length+1];

  for (int i=0; i<length; i++)
  {
    c[i] = str[i];
  }

  a.size = length;
  a.p = c;
  return a;
}

void print_find(Index* p_seek, Nucleotide* n_compare, int i, int k, int nucl_read, Nucleotide* n_seek, const bool& analysis)
{
  //p_seek is the pointer of the Index being modified
  //n_compare is the pointer to the last element modified
  //i is the position of the region of interest
  //k is the length of the region of interest minus 1
  //nucl_read is the total number of nucleotides in sequence
  //n_seek pointer to the first element of interest

  Nucleotide *last_10_pointer, *n_show;
  if (analysis)
    cout << "Sequence found in:\t" << p_seek->get_filename() << endl;

  cout << "base pair positions:\t[" << i << ":" << i+k  << "]"<< endl;
  cout << "base pair length:\t" << k+1 << endl;

  last_10_pointer = p_seek->get_n_head();

  //in case the position of the region of interest is less than ten, the previous characters are not 10, but less
  if (i<10)
  {
    cout << "prev " << i << " base pairs:\t";
    for(int j=0; j<i; j++)
    {
      cout << last_10_pointer->getN();
      last_10_pointer=last_10_pointer->getNextN();
    }
    cout << endl;
  }
  else
  {
    //if not, get pointers to the first element of the 10 previous nucleotides
    for(int j=0; j<(i-10); j++)
      last_10_pointer=last_10_pointer->getNextN();

    cout << "prev 10 base pairs:\t";
    for(int j=0; j<10; j++)
    {
      cout << last_10_pointer->getN();
      last_10_pointer=last_10_pointer->getNextN();
    }
    cout << endl;
  }

  cout << "region of interest:\t";
  for(int i=0; i<k+1; i++)
  {
    cout << n_seek->getN();
    n_seek = n_seek->getNextN();
  }

  n_show = n_compare->getNextN();

  //check if modified region is less that 10 nucleotides away from the end of the sequence.
  //If it is, less than 10 characters will be printed.
  if(nucl_read-k-i-2<10)
  {
    cout << "\nnext " << nucl_read-k-i-1<< " base pairs:\t";
    for(int j=0; j<nucl_read-k-i-1; j++)
    {
      cout << n_show->getN();
      n_show=n_show->getNextN();
    }
    cout << endl;
  }
  else
  {
    cout << "\nnext 10 base pairs:\t";
    for(int j=0; j<10; j++)
    {
      cout << n_show->getN();
      n_show=n_show->getNextN();
    }
    cout << endl;
  }
}

void option2_2(const string& filename, DNADatabase& dna_db)
{
  string user_input;
  bool good_file;
  bool show=true;
  TicToc time;

  while(show)
  {
    cout << "\nFind DNA sequence by file" << endl;
    cout << "\nEnter FASTA file to search (eg, a.fna):" << endl;
    getline(cin, user_input);

    time.tic();

    good_file = check_file(user_input);

    if (good_file)
    {
      dna_db.find_from_file(user_input, filename, false, 0);

      time.toc();
      cout << endl;
      cout << time << endl;

      show = false;
    }
    else
    {
      cout << "\nPlease try again" << endl;
    }
  }
}



void option2_3(const string& filename, DNADatabase& dna_db, bool& changes_done)
{
  string_for_class sequence;
  string user_input, str_pos;
  bool show=true;
  bool good_number;
  int position;
  TicToc time;

  cout << "\nAdd DNA sequence by input" << endl;
  cout << "\nEnter the DNA sequence to add:" << endl;
  getline(cin, user_input);
  while(show)
  {
    cout << "\nEnter a base pair position:" << endl;
    getline(cin, str_pos);
    time.tic();
    good_number = check_string_is_num(str_pos);
    if(good_number)
    {
      //get number from string
      position = stoi(str_pos);
      if (position<0)
        cout << "\nThe number selected is out of range. Please try again." << endl;
      else
        show = false;
    }
    else
      cout << "\nRemember to enter numbers only." << endl;
  }

  //get character array from string
  sequence = string_to_chars(user_input);

  dna_db.add_NUCL(sequence, position, filename, changes_done);
  time.toc();
  cout << endl;
  cout << time << endl;

  delete[] sequence.p;
}

void option2_4(const string& filename, DNADatabase& dna_db, bool& changes_done)
{
  string_for_class sequence;
  string user_input, str_pos;
  bool show_pos=true;
  bool show_file=true;
  bool good_file=false;;
  bool good_number;
  int position;
  TicToc time;

  cout << "\nAdd DNA sequence by file" << endl;
  while(show_file)
  {
    cout << "\nEnter file to add:" << endl;
    getline(cin, user_input);
    good_file = check_file(user_input);
    if (good_file)
      show_file =false;
    else cout << "\nPlease try again." << endl;
  }


  while(show_pos)
  {
    cout << "\nEnter a base pair position:" << endl;
    getline(cin, str_pos);
    time.tic();
    good_number = check_string_is_num(str_pos);
    if(good_number)
    {
      //get number from string
      position = stoi(str_pos);
      if (position<0)
        cout << "\nThe number selected is out of range. Please try again." << endl;
      else
        show_pos = false;
    }
    else
      cout << "\nRemember to enter numbers only." << endl;
  }

  dna_db.add_NUCL_from_file(user_input, position, filename, changes_done);

  time.toc();
  cout << endl;
  cout << time << endl;

  delete[] sequence.p;
}

void option2_5(const string& filename, DNADatabase& dna_db, bool& changes_done)
{
  string str_pos, str_length;
  bool show=true;
  bool good_number;
  int position;
  int length;
  TicToc time;

  cout << "\nDelete DNA sequence by input" << endl;

  while(show)
  {
    cout << "\nEnter a base pair position:" << endl;
    getline(cin, str_pos);
    good_number = check_string_is_num(str_pos);
    if(good_number)
    {
      //get number from string
      position = stoi(str_pos);
      if (position<0)
        cout << "\nThe number selected is out of range. Please try again." << endl;
      else
        show = false;
    }
    else
      cout << "\nRemember to enter numbers only." << endl;
  }

  show = true;

  while(show)
  {
    cout << "\nEnter a base pair length:" << endl;
    getline(cin, str_length);
    time.tic();
    good_number = check_string_is_num(str_length);
    if(good_number)
    {
      //get number from string
      length = stoi(str_length);
      if (length<0)
        cout << "\nThe number selected is out of range. Please try again." << endl;
      else
        show = false;
    }
    else
      cout << "\nRemember to enter numbers only." << endl;
  }

  dna_db.deleteN(position, length, filename, changes_done);

  time.toc();

  cout << endl;
  cout << time << endl;

}

void option2_6(const string& filename, DNADatabase& dna_db, bool& changes_done)
{
  string str_pos, str_length, str_seq;
  bool show=true;
  bool good_number, good_file;
  int position;
  int length;
  TicToc time;
  string_for_class sequence;

  cout << "\nReplace DNA sequence by input" << endl;

  cout << "\nEnter the DNA sequence to add" << endl;
  getline(cin, str_seq);

  //create array from string
  sequence = string_to_chars(str_seq);

  while(show)
  {
    cout << "\nEnter a base pair position:" << endl;
    getline(cin, str_pos);
    good_number = check_string_is_num(str_pos);
    if(good_number)
    {
      //get number from string
      position = stoi(str_pos);
      if (position<0)
        cout << "\nThe number selected is out of range. Please try again." << endl;
      else
        show = false;
    }
    else
      cout << "\nRemember to enter numbers only." << endl;
  }

  show = true;

  while(show)
  {
    cout << "\nEnter a base pair length:" << endl;
    getline(cin, str_length);
    time.tic();
    good_number = check_string_is_num(str_length);
    if(good_number)
    {
      //get number from string
      length = stoi(str_length);
      if (length<0)
        cout << "\nThe number selected is out of range. Please try again." << endl;
      else
        show = false;
    }
    else
      cout << "\nRemember to enter numbers only." << endl;
  }

  dna_db.change(position, length, sequence, filename, changes_done);
  time.toc();

  cout << endl;
  cout << time << endl;

  delete[] sequence.p;
}

void option2_7(const string& filename, DNADatabase& dna_db, bool& changes_done)
{
  string str_pos, str_length, file_to_look;
  bool show=true;
  bool good_number;
  bool good_file;
  bool show_file=true;
  int position;
  int length;
  TicToc time;

  cout << "\nReplace DNA sequence by file" << endl;
  while (show_file)
  {
    cout << "\nEnter file to add" << endl;
    getline(cin, file_to_look);
    good_file = check_file(file_to_look);
    if (!good_file)
      cout << "\nPlease try again." << endl;
    else
      show_file = false;
  }

  while(show)
  {
    cout << "\nEnter a base pair position:" << endl;
    getline(cin, str_pos);
    good_number = check_string_is_num(str_pos);
    if(good_number)
    {
      //get number from string
      position = stoi(str_pos);
      if (position<0)
        cout << "\nThe number selected is out of range. Please try again." << endl;
      else
        show = false;
    }
    else
      cout << "\nRemember to enter numbers only." << endl;
  }

  show = true;

  while(show)
  {
    cout << "\nEnter a base pair length:" << endl;
    getline(cin, str_length);
    time.tic();
    good_number = check_string_is_num(str_length);
    if(good_number)
    {
      //get number from string
      length = stoi(str_length);
      if (length<0)
        cout << "\nThe number selected is out of range. Please try again." << endl;
      else
        show = false;
    }
    else
      cout << "\nRemember to enter numbers only." << endl;
  }

  dna_db.change_from_file(position, length, file_to_look, filename, changes_done);
  time.toc();

  cout << endl;
  cout << time << endl;
}


void option2_8(DNADatabase& dna_db, string filename, const bool& changes_done)
{
  string user_input;
  TicToc time;
  cout << "\nSave edited DNA sequence" << endl;

  if (!changes_done)
    cout << "\nYou have not edited any sequence yet." << endl;
  else
  {
    cout << "\nEnter desired filename for your edited sequence (eg, b.fna)" << endl;
    getline(cin, user_input);
    time.tic();
    dna_db.save(filename, user_input);
    time.toc();
    cout << endl;
    cout << time << endl;
  }
}

void print_before(Nucleotide* n_seek, int position)
{
  //n_seek is the pointer to the first element of the sequence being modified

  char nucl;

  //check if the position of the region of interest if less than 10 characters away from the start of the sequence.
  //If it is, less than 10 characters will be printed.
  if (position<10)
  {
    cout << "prev " << position << " base pairs:\t";
    for(int i=0; i<position; i++)
    {
      nucl = n_seek->getN();
      cout << nucl;
      n_seek = n_seek->getNextN();
    }
    cout << endl;
  }
  else
  {
    cout << "prev 10 base pairs:\t";
    for(int i=0; i<position-10; i++)
        n_seek = n_seek->getNextN();

    for(int j=0; j<10; j++)
    {
      nucl = n_seek->getN();
      cout << nucl;
      n_seek = n_seek->getNextN();
    }
    cout << endl;
  }
}

void print_after(Nucleotide* n_seek, int position_last, int nucl_number)
{
  //n_seek points to the character after the region of interest

  //If the end of the region of interest is less than 10 characters away, less than 10 characters will be printed
  if (nucl_number-position_last<10)
  {
    cout << "next " << nucl_number-position_last << " base pairs:\t";
    for (int i=0; i<nucl_number-position_last; i++)
      {
        cout << n_seek->getN();
        n_seek = n_seek->getNextN();
      }
  }
  else
  {
    cout << "next 10 base pairs:\t";
    for(int i=0; i<10; i++)
    {
      cout << n_seek->getN();
      n_seek = n_seek->getNextN();
    }
  }
  cout << endl;
}

void option3(vector<string>& filename_list, DNADatabase& dna_db)
{
  string user_input;
  TicToc time;
  bool show_file = true;
  bool good_file;

  if (filename_list.empty())
    cout << "\nYou have not added any sequences to the database." << endl;
  else
  {
    cout << "\nAnalyse the DNA database" << endl;
    while(show_file)
    {
      cout << "\nEnter file name with .fna extension." << endl;
      getline(cin, user_input);
      good_file = check_file(user_input);
      if (good_file)
        show_file = false;
      else
        cout << "\nPlease try again." << endl;
    }
      time.tic();
      dna_db.analysis(user_input, filename_list);
      time.toc();
      cout << endl;
      cout << time << endl;
  }
}
