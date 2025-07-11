@def hassim=true;
@def title = "Chen solver"
@def t0 = "0"
@def tf = "120"
@def var1 = "x"
@def var1_value = "-0.1"
@def var2 = "y"
@def var2_value = "0.5"
@def var3 = "z"
@def var3_value = "-0.6"
@def params = ["a", "b", "c"]
@def a = "2"
@def b = "1"
@def c = "-0.5"

~~~
<head>
<script src="/libs/common/generateTableXYZ.js"></script>
</head>
~~~

~~~
        <div>
            This webpage uses the <a href='/RKF45/' link='_blank'>Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)</a> to approximate the solution to <a href='https://en.wikipedia.org/wiki/Multiscroll_attractor' link='_blank'>Chen equations</a>:

            \[
            \begin{aligned}
            \dfrac{dx}{dt} &= a (y-x) \\
            \dfrac{dy}{dt} &= x(c-a-z) + cy \\
            \dfrac{dz}{dt} &= xy - b z.
            \end{aligned}
            \]

        </div>

        <!--A form for users to enter in all the parameters of the problem-->
        <form name="requiredData">
            <table>
                <tr style="border: 0px solid black; padding: 0px;">
                    <th style="border: 1px solid black; padding: 0px;">Parameter</th>
                    <th style="border: 1px solid black; padding: 0px;">Value</th>
                    <th style="border: 1px solid black; padding: 0px;">Explanation</th>
                </tr>
                <tr>
                    <td><label for="a">\(a\):</label></td>
                    <td><input type="Number"
                        id="a"
                        name="a"
                        value="40"></td>
                    <td>Problem parameter.</td>
                </tr>
                <tr>
                    <td><label for="b">\(b\):</label></td>
                    <td><input type="Number"
                        id="b"
                        name="b"
                        value="3"></td>
                    <td>Problem parameter.</td>
                </tr>
                <tr>
                    <td><label for="c">\(c\):</label></td>
                    <td><input type="Number"
                        id="c"
                        name="c"
                        value="28"></td>
                    <td>Problem parameter.</td>
                </tr>
                {{ render_form(params) }}
                {{ insert attractor_init_coords.html }}
                {{ insert attractor_data.html }}
            </table>
        </form>

        <!--Buttons-->
        {{ insert attractor_button_table.html }}

        <!--Where the table and plot goes-->
        {{ insert attractor_output.html }}
~~~