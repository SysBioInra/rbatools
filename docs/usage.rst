Usage
=====

.. _installation:

Installation
------------

To use rbatools, first install it using pip:

.. code-block:: console

    pip install rbatools

Usage remark
------------

To use rbatools, import it in Python

.. code-block:: console

    import rbatools


Learning materials and further information
------------------------------------------
For further information, please visit `our website <https://rba.inrae.fr>` and for tutorials our `GitHub <https://github.com/SysBioInra/rbatools>`.

Usage remark
------------

rbatools is based on linear-programming and therefore requires a linear-optimization package.
We strongly suggest using the (propriatary) CPLEX tool, which is available under a free academic license.

Nonetheless we included functionality with the freely available swiglpk solver, however this sometimes is not able to obtain satisfactory solutions.

