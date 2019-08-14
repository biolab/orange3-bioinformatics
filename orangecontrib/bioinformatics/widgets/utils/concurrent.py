""" Multithreading widgets with QThreadPool """

from AnyQt.QtCore import QObject, QRunnable, pyqtSlot, pyqtSignal


class Signals(QObject):
    progress = pyqtSignal()
    partial_result = pyqtSignal(object)
    result = pyqtSignal(object)
    finished = pyqtSignal()
    error = pyqtSignal(Exception)


class Worker(QRunnable):
    """ Worker runnable
    """

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = Signals()

        is_callback = kwargs.get("progress_callback", None)
        is_partial = kwargs.get("partial_result", None)

        if is_callback:
            self.kwargs['progress_callback'] = self.signals.progress

        if is_partial:
            self.kwargs['partial_result'] = self.signals.partial_result

    @pyqtSlot()
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception as e:
            self.signals.error.emit(e)
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()
