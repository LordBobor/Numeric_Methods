import org.math.plot.Plot2DPanel;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class Parabolic_Graphics {
    private  JTextField a_textField;
    private  JTextField b_textField;
    private JTextField c_textField;
    private JFrame frame;
    private JPanel Test;
    private JButton button1;
    private JCheckBox точноеРешениеCheckBox;
    private JCheckBox явныйМетодCheckBox;
    private JCheckBox неявныйМетодCheckBox;
    private Parabolic solver;
    Plot2DPanel plot;

    public static void main(String args[]){
        new Parabolic_Graphics();
    }

    public Parabolic_Graphics() {
        initSolver();
        initFrames();
        initListeners();
    }

    private void initFrames() {


        plot = new Plot2DPanel();


        JPanel panel = new JPanel();
        panel.setLayout(new BorderLayout());
        panel.add(plot);

        JPanel container = new JPanel();
        container.setLayout(new BoxLayout(container, BoxLayout.X_AXIS));

        container.add(Test);
        container.add(panel);



        frame = new JFrame("Main");
        frame.add(container);
        // Init frame
       // frame.setContentPane(Test);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);
        frame.setSize(new Dimension(800, 600));


        button1.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                plot.removeAllPlots();
                solver.setA(a_textField.getText());
                solver.setB(b_textField.getText());
                solver.setC(c_textField.getText());

                double Xn = Math.PI / 2;
                double h = Xn / 20;

                double[] x_0 = new double[21];
                double[] y_0 = new double[21];

                for (int i = 0; i <= 20; ++i) {
                    x_0[i] = i * h;
                    y_0[i] = solver.exactSolution(x_0[i], 1);
                }

                solver.Solve(0);
                double[] x_2 = solver.getX();
                double[] y_2 = solver.getRes();

                plot.addLinePlot("exact plot", Color.red, x_0, y_0);
                plot.addLinePlot("implicit plot", Color.blue, x_2, y_2);

            }
        });

    }

    private void initSolver() {
        solver = new Parabolic(0.02, 0, -0.04, 20, 200);
    }

    private void initListeners() {
        a_textField.getDocument().addDocumentListener(new myDocumentListener());
        b_textField.getDocument().addDocumentListener(new myDocumentListener());
        c_textField.getDocument().addDocumentListener(new myDocumentListener());
    }

    private void createUIComponents() {
        // TODO: place custom component creation code here
    }


    public class myDocumentListener implements DocumentListener {
        @Override
        public void insertUpdate(DocumentEvent e) {
           action();
        }
        @Override
        public void removeUpdate(DocumentEvent e) {
            action();
        }
        @Override
        public void changedUpdate(DocumentEvent e) {
            action();
        }
        public void action(){
            //a_textField.setText(a_textField.getText());
            //solver.setN("10");
            //solver.setK("10");
        }
    }
}
