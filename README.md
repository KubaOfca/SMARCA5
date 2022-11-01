<div style="text-align: center; font-size: 30px">SMARCA5</div>

---

# Description

Main goal of this application is based on making an alignment of any sequence (like protein, DNA, RNA ect.) 
to its subsequence (f.e. peptides). Visualisation of this proces is done by parsing an output into HTML and CSS,
so you only need browser to see the result in clear way. In addition, ```SMARCA5``` calculates basic statistic of 
any subsequence which are displayed on mouse hoover. If you have result from more than one experiment you
will be able to switch a view between them (you can even divide it on group and samples).

---

# Requirements

## File formats

### Sequence

The sequence data should be ```.fasta``` or simple ```.txt``` file format.

### Subsequence

The best way to represent subsequence data is ```.csv``` or ```.xlsx``` file format.
Column name are specific and must look like this:

| Sequence | Experiment | Protein            |
|----------|------------|--------------------|
| MMAVA    | Exp1       | SMARCA5            |
| KKKMV    | Exp2       | SMARCA5; SMARCA5.1 |
| ........ | .......... | .................. |


In ```Sequence``` column you should put your subsequence that you obtain from experiment

In ```Experiment``` column you should put your experiment label in given format (like below):

```<groupName>_<treatmentType><sampleNumber>```

In ```Protein``` column you should put information about form which protein this subsequence come (it can be found 
in databases like ```UnitProt```). If subsequence is not unique (appears in more than one protein) put ```;```(semicolon)
delimiter between them.

---

# Display options

```SMARCA5``` provides two style format to change:
1. Line length
   - it allows you to display more aminoacids in one line
2. Interline (given in pixels)
   - it allows you to increase or decrease space between lines

# DEMO

1. Load ```.fasta``` and ```.csv``` file.
2. Adjust ```line lenght``` and ```interline``` options.

![image](https://user-images.githubusercontent.com/61982713/194770278-7e5d5598-e554-4e24-a1ba-19fd20585bf9.png)

3. Click ```Anlayze``` button.
4. After this new browser tab should pop up.

![image](https://user-images.githubusercontent.com/61982713/194770292-a15347f1-6840-40b0-83de-df5e0dd67afe.png)

5. You can display subseqence details by mouse hover

![image](https://user-images.githubusercontent.com/61982713/194770332-68250e07-4c8a-4c9a-acd9-1dca0d7dffd4.png)

6. Menu (on right side) allows you to switch a view between different experiments

(Photo)


# Installation

## Clone repo

1. Clone git repo
```bash
git clone https://github.com/KubaOfca/SMARCA5.git
git cd SMARCA5
```
2. [Install Poetry](https://python-poetry.org/docs/#installation) (python packaging and dependency management)
3. Use Poetry to create and activate the environment by typing:
```bash
poetry shell
```
4. Install `SMARCA5` dependencies by typing:
```bash
poetry install
```

## Download .exe version
