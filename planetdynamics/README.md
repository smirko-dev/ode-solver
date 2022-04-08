# Planet dynamics

## File format

```sh
10
```
Number of planets (uint32_t)

```sh
100.0 100.0 1.0 1.0 5.0 0.1
```
Planet data containing the name (string) position [x] (float_t), velocity [dx] (float_t), mass (float_t) and radius (float_t)

## Usage

```sh
pd planets.dat
```