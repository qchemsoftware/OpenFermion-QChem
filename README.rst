Q-CHEM interfaced with OpenFermion
=================

`OpenFermion <http://openfermion.org>`__ is an open source library (licensed under Apache 2) for compiling and analyzing quantum algorithms which simulate fermionic systems.
This plugin library allows the electronic structure package `Q-CHEM <https://www.q-chem.com>`__ to interface with OpenFermion.
Users who want to use OpenFermion/Q-CHEM interface need to obtain Q-CHEM software and agree to the Q-CHEM lisence separately.

Requirements
------------
- Q-CHEM 5.4+
- OpenFermion 1.3+
- all requirements for the aforementioned programs

Installation
------------
To use OpenFermion/Q-CHEM interface:

- Install Q-CHEM: `https://www.q-chem.com/install <https://www.q-chem.com/install>`__
- Install OpenFermion: `https://github.com/quantumlib/OpenFermion <https://github.com/quantumlib/OpenFermion>`__
- Install OpenFermion/Q-CHEM:

.. code-block:: bash

  git clone https://github.com/qchemsoftware/OpenFermion-QChem.git
  cd OpenFermion-QChem
  python -m pip install -e .

PyPI releases soon to be available.

Please take a look at the `ipython notebook demo <https://github.com/qchemsoftware/OpenFermion-QChem/blob/main/examples/Openfermion-Qchem.ipynb>`__ and general description of `OpenFermion/Q-CHEM v1.0 <https://github.com/qchemsoftware/OpenFermion-QChem/blob/main/examples/OpenFermion-QCHEM_v1.0.pdf>`__.

Authors
-------

`Anna I. Krylov <https://iopenshell.usc.edu/>`__ (University of Southern California),
`Yongbin Kim <https://github.com/yongbinkim-chemist>`__ (University of Southern California),
`Evgeny Epifanovsky <https://www.q-chem.com/>`__ (Q-CHEM)

How to cite
-----------
When using OpenFermion/Q-CHEM for research projects, please cite:

    Jarrod R. McClean, Kevin J. Sung, Ian D. Kivlichan, Yudong Cao,
    Chengyu Dai, E. Schuyler Fried, Craig Gidney, Brendan Gimby,
    Pranav Gokhale, Thomas Häner, Tarini Hardikar, Vojtĕch Havlíček,
    Oscar Higgott, Cupjin Huang, Josh Izaac, Zhang Jiang, Xinle Liu,
    Sam McArdle, Matthew Neeley, Thomas O'Brien, Bryan O'Gorman, Isil Ozfidan,
    Maxwell D. Radin, Jhonathan Romero, Nicholas Rubin, Nicolas P. D. Sawaya,
    Kanav Setia, Sukin Sim, Damian S. Steiger, Mark Steudtner, Qiming Sun,
    Wei Sun, Daochen Wang, Fang Zhang and Ryan Babbush.
    *OpenFermion: The Electronic Structure Package for Quantum Computers*.
    `arXiv:1710.07629 <https://arxiv.org/abs/1710.07629>`__. 2017.

as well as

    Evgeny Epifanovsky, Andrew T. B. Gilbert, Xintian Feng, Joonho Lee, Yuezhi Mao,
    Narbe Mardirossian, Pavel Pokhilko, Alec F. White, Marc P. Coons, Adrian L. Dempwolff,
    Zhengting Gan, Diptarka Hait, Paul R. Horn, Leif D. Jacobson, Ilya Kaliman, Jörg Kussmann,
    Adrian W. Lange, Ka Un Lao, Daniel S. Levine, Jie Liu, Simon C. McKenzie, Adrian F. Morrison,
    Kaushik D. Nanda, Felix Plasser, Dirk R. Rehn, Marta L. Vidal, Zhi-Qiang You, Ying Zhu,
    Bushra Alam, Benjamin J. Albrecht, Abdulrahman Aldossary, Ethan Alguire, Josefine H. Andersen,
    Vishikh Athavale, Dennis Barton, Khadiza Begam, Andrew Behn, Nicole Bellonzi, Yves A. Bernard,
    Eric J. Berquist, Hugh G. A. Burton, Abel Carreras, Kevin Carter-Fenk, Romit Chakraborty,
    Alan D. Chien, Kristina D. Closser, Vale Cofer-Shabica, Saswata Dasgupta, Marc de Wergifosse,
    Jia Deng, Michael Diedenhofen, Hainam Do, Sebastian Ehlert, Po-Tung Fang, Shervin Fatehi,
    Qingguo Feng, Triet Friedhoff, James Gayvert, Qinghui Ge, Gergely Gidofalvi, Matthew Goldey,
    Joe Gomes, Cristina E. González-Espinoza, Sahil Gulania, Anastasia O. Gunina,
    Magnus W. D. Hanson-Heine, Phillip H. P. Harbach, Andreas Hauser, Michael F. Herbst,
    Mario Hernández Vera, Manuel Hodecker, Zachary C. Holden, Shannon Houck, Xunkun Huang,
    Kerwin Hui, Bang C. Huynh, Maxim Ivanov, Ádám Jász, Hyunjun Ji, Hanjie Jiang, Benjamin Kaduk,
    Sven Kähler, Kirill Khistyaev, Jaehoon Kim, Gergely Kis, Phil Klunzinger, Zsuzsanna Koczor-Benda,
    Joong Hoon Koh, Dimitri Kosenkov, Laura Koulias, Tim Kowalczyk, Caroline M. Krauter, Karl Kue,
    Alexander Kunitsa, Thomas Kus, István Ladjánszki, Arie Landau, Keith V. Lawler, Daniel Lefrancois,
    Susi Lehtola, Run R. Li, Yi-Pei Li, Jiashu Liang, Marcus Liebenthal, Hung-Hsuan Lin,
    You-Sheng Lin, Fenglai Liu, Kuan-Yu Liu, Matthias Loipersberger, Arne Luenser, Aaditya Manjanath,
    Prashant Manohar, Erum Mansoor, Sam F. Manzer, Shan-Ping Mao, Aleksandr V. Marenich,
    Thomas Markovich, Stephen Mason, Simon A. Maurer, Peter F. McLaughlin, Maximilian F. S. J. Menger,
    Jan-Michael Mewes, Stefanie A. Mewes, Pierpaolo Morgante, J. Wayne Mullinax,
    Katherine J. Oosterbaan, Garrette Paran, Alexander C. Paul, Suranjan K. Paul, Fabijan Pavošević,
    Zheng Pei, Stefan Prager, Emil I. Proynov, Ádám Rák, Eloy Ramos-Cordoba, Bhaskar Rana,
    Alan E. Rask, Adam Rettig, Ryan M. Richard, Fazle Rob, Elliot Rossomme, Tarek Scheele,
    Maximilian Scheurer, Matthias Schneider, Nickolai Sergueev, Shaama M. Sharada,
    Wojciech Skomorowski, David W. Small, Christopher J. Stein, Yu-Chuan Su, Eric J. Sundstrom,
    Zhen Tao, Jonathan Thirman, Gábor J. Tornai, Takashi Tsuchimochi, Norm M. Tubman,
    Srimukh Prasad Veccham, Oleg Vydrov, Jan Wenzel, Jon Witte, Atsushi Yamada, Kun Yao, Sina Yeganeh,
    Shane R. Yost, Alexander Zech, Igor Ying Zhang, Xing Zhang, Yu Zhang, Dmitry Zuev,
    Alán Aspuru-Guzik, Alexis T. Bell, Nicholas A. Besley, Ksenia B. Bravaya, Bernard R. Brooks,
    David Casanova, Jeng-Da Chai, Sonia Coriani, Christopher J. Cramer, György Cserey,
    A. Eugene DePrince III, Robert A. DiStasio Jr., Andreas Dreuw, Barry D. Dunietz,
    Thomas R. Furlani, William A. Goddard III, Sharon Hammes-Schiffer, Teresa Head-Gordon,
    Warren J. Hehre, Chao-Ping Hsu, Thomas-C. Jagau, Yousung Jung, Andreas Klamt, Jing Kong,
    Daniel S. Lambrecht, WanZhen Liang, Nicholas J. Mayhall, C. William McCurdy, Jeffrey B. Neaton,
    Christian Ochsenfeld, John A. Parkhill, Roberto Peverati, Vitaly A. Rassolov, Yihan Shao,
    Lyudmila V. Slipchenko, Tim Stauch, Ryan P. Steele, Joseph E. Subotnik, Alex J. W. Thom,
    Alexandre Tkatchenko, Donald G. Truhlar, Troy Van Voorhis, Tomasz A. Wesolowski,
    K. Birgitta Whaley, H. Lee Woodcock III, Paul M. Zimmerman, Shirin Faraji, Peter M. W. Gill,
    Martin Head-Gordon, John M. Herbert, and Anna I. Krylov.
    *Software for the frontiers of quantum chemistry: An overview of developments in the Q-Chem 5 package*.
    `DOI: 10.1063/5.0055522 <https://aip.scitation.org/doi/10.1063/5.0055522>`__.
    2017.

We are happy to include future contributors as authors on later OpenFermion/Q-CHEM releases.
