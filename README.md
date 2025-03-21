# Control

Some code snippets related to control theory, transfer functions, and such. At the moment, there is only one Matlab function to discretize continuous dynamic models (transfer functions or state-space models), `cont2disc`. The available discretization methods are

- `backward` - backward difference
- `forward`  - forward difference
- `euler`    - same as `forward`
- `tustin`   - tustin method
- `bilinear` - same as `tustin`
- `impulse`  - impulse-invariance
- `step`     - step-invariance
- `ramp`     - ramp-invariance
- `zoh`      - zero-order hold

Two other functions, needed by `cont2disc`, are also included: `tf2sym` and `sym2tf`.